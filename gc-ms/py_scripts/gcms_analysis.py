#!/usr/bin/env python3
"""
Data analysis for GC/MS QC chromatograms: from an mzML file to matched target
compounds and their quality metrics.

This module is self-contained (no plotting, no report writing) and can be used
on its own:

    from gcms_analysis import load_expected_compounds, analyze_mzml

    compounds = load_expected_compounds('compounds.csv')   # name -> (RT in seconds, external ID)
    analysis = analyze_mzml('sample.mzML', compounds)
    print(analysis.metrics['n_found'], analysis.metrics['sn_quantiles'])

Pipeline (see analyze_chromatogram):
  1. Compute the TIC of an MSExperiment as an MSChromatogram
  2. Pick chromatographic peaks (chromatogram_peak_picker)
  3. Match the picked peaks against the expected compounds (name -> RT), with
     RT-bias correction
  4. Summarize the matches into run-level metrics

It also owns the compound list file format (see load_expected_compounds) and the
per-file TSV exports of its results.
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path
from typing import Any, NamedTuple, TypedDict

from pyopenms import (
    MSChromatogram,
    MSExperiment,
    MSSpectrum,
    MzMLFile,
    Peak1D,
    SpectrumAlignment,
)

from chromatogram_peak_picker import pick_chromatographic_peaks

class PeakInfo(TypedDict):
    """One picked chromatographic peak, as returned by pick_chromatographic_peaks."""
    rt: float
    intensity: float
    fwhm: float | None
    sn: float

class ExpectedCompound(NamedTuple):
    """One expected (target) compound of the compound list, which is keyed by name."""
    rt: float           # expected retention time, in seconds
    external_id: str    # identifier in an external database, e.g. 'CHEBI:15882'

class MatchResult(TypedDict):
    """One expected compound matched (or not) against a picked peak."""
    compound_name: str
    external_id: str
    expected_rt: float
    observed_rt: float | None
    delta_rt: float | None
    intensity: float | None
    fwhm: float | None
    sn: float | None

class RunMetrics(TypedDict):
    """
    Run-level summary of one file's matched compounds: the quartiles (Q1, Q2, Q3)
    of delta RT / S/N / FWHM -- Q2 being the median and Q3 - Q1 the IQR -- the mean
    intensity, and how many of the expected compounds were found.
    """
    delta_rt_quantiles: tuple[float, float, float] | None
    sn_quantiles: tuple[float, float, float] | None
    fwhm_quantiles: tuple[float, float, float] | None
    mean_intensity: float | None
    n_found: int
    n_expected: int

class FileAnalysis(NamedTuple):
    """Everything analyze_chromatogram derived from one chromatogram."""
    chromatogram: MSChromatogram          # raw TIC
    smoothed_chromatogram: MSChromatogram # smoothed TIC the peaks were picked on
    noise_levels: list[float]             # noise estimate per point of the smoothed TIC
    observed_peaks: list[PeakInfo]        # all picked peaks
    matches: list[MatchResult]            # one entry per expected compound ([] if no peak was picked)
    metrics: RunMetrics                   # run-level summary of `matches`

# Column layout of the compound list CSV: a name column, an external database
# identifier column (required by the mzQC per-compound tables, see gcms_mzqc) plus
# exactly one RT column, which fixes the unit. Comment (#) and blank lines are ignored.
COMPOUND_NAME_COLUMN = 'compound_name'
COMPOUND_EXTERNAL_ID_COLUMN = 'external_id'
COMPOUND_RT_COLUMNS: dict[str, float] = {  # column name -> factor converting it to seconds
    'rt_seconds': 1.0,
    'rt_minutes': 60.0,
}
COMPOUNDS_TEMPLATE = """\
# Expected (target) compounds for GC/MS QC, one per line.
#
# Required columns: '{name_col}', '{id_col}' and exactly one of '{rt_cols}',
# which determines the retention time unit. '{id_col}' is the compound's identifier
# in an external database, e.g. a ChEBI accession. Lines starting with '#' are ignored.
# Replace the two example compounds below with your own.
{name_col},{id_col},rt_minutes
n-Nonane,CHEBI:32892,3.103
Phenol,CHEBI:15882,3.5
"""

def write_compounds_template(output_filename: str | Path) -> None:
    """
    Write a commented compound list template with two example compounds, meant to
    be adapted by the user and read back by load_expected_compounds.
    """
    with open(output_filename, 'w', newline='') as f:
        f.write(COMPOUNDS_TEMPLATE.format(name_col=COMPOUND_NAME_COLUMN,
                                          id_col=COMPOUND_EXTERNAL_ID_COLUMN,
                                          rt_cols="', '".join(COMPOUND_RT_COLUMNS)))

def load_expected_compounds(filename: str | Path) -> dict[str, ExpectedCompound]:
    """
    Read the expected (target) compounds from a CSV file.

    The file needs a header with a '<COMPOUND_NAME_COLUMN>' column, a
    '<COMPOUND_EXTERNAL_ID_COLUMN>' column (the compound's identifier in an external
    database, e.g. a ChEBI accession) and exactly one retention time column from
    COMPOUND_RT_COLUMNS ('rt_seconds' or 'rt_minutes'), which fixes the unit of all
    RTs. Blank lines and lines starting with '#' are ignored, so the file can carry
    comments (see write_compounds_template).

    Returns:
        Dict mapping compound name -> ExpectedCompound, with the expected RT in
        *seconds* (the unit used throughout this module, as in pyOpenMS).

    Raises:
        ValueError: on a missing/ambiguous header, an empty compound name or
        external ID, an unparsable or negative RT, a duplicate compound name, or
        an empty compound list.
    """
    filename = Path(filename)
    with open(filename, newline='') as f:
        rows = [line for line in f if line.strip() and not line.lstrip().startswith('#')]

    reader = csv.DictReader(rows)
    columns = reader.fieldnames or []
    for column in (COMPOUND_NAME_COLUMN, COMPOUND_EXTERNAL_ID_COLUMN):
        if column not in columns:
            raise ValueError(f"{filename}: missing required column '{column}' "
                             f"(found: {', '.join(columns) or 'nothing'})")
    rt_columns = [c for c in columns if c in COMPOUND_RT_COLUMNS]
    if len(rt_columns) != 1:
        raise ValueError(f"{filename}: expected exactly one retention time column "
                         f"({' or '.join(COMPOUND_RT_COLUMNS)}), found "
                         f"{len(rt_columns)}: {', '.join(rt_columns) or 'none'}")
    rt_column = rt_columns[0]
    to_seconds = COMPOUND_RT_COLUMNS[rt_column]

    compounds: dict[str, ExpectedCompound] = {}
    for line_no, row in enumerate(reader, start=2):
        name = (row[COMPOUND_NAME_COLUMN] or '').strip()
        external_id = (row[COMPOUND_EXTERNAL_ID_COLUMN] or '').strip()
        raw_rt = (row[rt_column] or '').strip()
        if not name:
            raise ValueError(f"{filename}, compound {line_no - 1}: empty compound name")
        if name in compounds:
            raise ValueError(f"{filename}, compound {line_no - 1}: duplicate compound name '{name}'")
        if not external_id:
            raise ValueError(f"{filename}, compound '{name}': empty {COMPOUND_EXTERNAL_ID_COLUMN}")
        try:
            rt = float(raw_rt)
        except ValueError:
            raise ValueError(f"{filename}, compound '{name}': "
                             f"'{raw_rt}' is not a valid {rt_column} value") from None
        if rt < 0:
            raise ValueError(f"{filename}, compound '{name}': negative retention time {rt}")
        compounds[name] = ExpectedCompound(rt=rt * to_seconds, external_id=external_id)

    if not compounds:
        raise ValueError(f"{filename}: no compounds found")
    return compounds

def read_mzml_file(filename: str | Path) -> MSExperiment:
    """Read an mzML file and return the MSExperiment."""
    exp = MSExperiment()
    MzMLFile().load(str(filename), exp)
    return exp

def compute_chromatogram(exp: MSExperiment) -> MSChromatogram:
    """
    Compute the chromatogram using MSExperiment::calculateTIC.
    Returns an MSChromatogram object directly.
    """
    return exp.calculateTIC()

def _rt_spectrum(rts: list[float]) -> MSSpectrum:
    """Build an MSSpectrum with the given (ascending) RTs stored as m/z, for use with SpectrumAlignment."""
    spectrum = MSSpectrum()
    for rt in rts:
        spectrum.push_back(Peak1D(rt, 0.0))  # intensity is irrelevant for alignment
    return spectrum

def match_peaks_to_expected(
    expected_compounds: dict[str, ExpectedCompound],
    observed_peaks: list[PeakInfo],
    rt_tolerance: float,
    rt_bias: float = 0.0,
) -> list[MatchResult]:
    """
    Match observed peaks to expected compounds (name + RT) by aligning their RT
    positions with pyOpenMS's SpectrumAlignment: both the expected and observed
    RTs are loaded as m/z values into synthetic MSSpectrum objects (sorted by RT),
    and SpectrumAlignment's banded alignment matches them within rt_tolerance,
    leaving spurious observed peaks and missing compounds unmatched.

    Args:
        expected_compounds: Dict mapping compound name -> ExpectedCompound (RT in seconds)
        rt_bias: seconds added to each expected RT before aligning (e.g. to correct
            a systematic RT drift). The reported 'delta_rt' stays relative to the
            true (unbiased) expected_rt.

    Returns:
        List of dictionaries with: 'compound_name', 'external_id', 'expected_rt', 'observed_rt', 'delta_rt', 'intensity', 'fwhm', 'sn'
    """
    names_and_compounds = sorted(expected_compounds.items(), key=lambda kv: kv[1].rt)
    observed_peaks = sorted(observed_peaks, key=lambda p: p['rt'])

    aligner = SpectrumAlignment()
    params = aligner.getParameters()
    params.setValue("tolerance", rt_tolerance)
    params.setValue("is_relative_tolerance", "false")
    aligner.setParameters(params)

    theoretical_spectrum = _rt_spectrum([compound.rt + rt_bias for _, compound in names_and_compounds])
    observed_spectrum = _rt_spectrum([peak['rt'] for peak in observed_peaks])
    matched_obs_idx = dict(aligner.getSpectrumAlignment(theoretical_spectrum, observed_spectrum))

    results: list[MatchResult] = []
    for i, (name, compound) in enumerate(names_and_compounds):
        peak = observed_peaks[matched_obs_idx[i]] if i in matched_obs_idx else None
        results.append({
            'compound_name': name,
            'external_id': compound.external_id,
            'expected_rt': compound.rt,
            'observed_rt': peak['rt'] if peak else None,
            'delta_rt': (peak['rt'] - compound.rt) if peak else None,
            'intensity': peak['intensity'] if peak else None,
            'fwhm': peak['fwhm'] if peak else None,
            'sn': peak['sn'] if peak else None,
        })

    return results

def match_with_rt_bias_correction(
    expected_compounds: dict[str, ExpectedCompound],
    observed_peaks: list[PeakInfo],
    rt_tolerance_bias_estimation: float = 6.0,
    rt_tolerance: float = 3.0,
) -> list[MatchResult]:
    """
    Match observed peaks to expected compounds, correcting for a systematic RT
    drift in two steps:
      1. Estimate the drift by aligning all expected compounds (uncorrected) against
         only the n most intense picked peaks overall (n = number of expected
         compounds, assumed to be exactly the n target compounds), using the looser
         rt_tolerance_bias_estimation, and taking the median signed delta RT of
         whatever matched as the RT bias.
      2. Align all picked peaks against the bias-corrected expected RTs for the
         final result, using the tighter rt_tolerance.

    rt_tolerance_bias_estimation must be looser than rt_tolerance: the bias can only
    ever be estimated from compounds whose (uncorrected) delta RT already falls
    within it, so a drift as large as or larger than rt_tolerance would otherwise
    never be detected in step 1 -- making step 2's correction a no-op exactly when
    it's needed most.
    """
    n = len(expected_compounds)
    top_peaks = sorted(observed_peaks, key=lambda p: p['intensity'], reverse=True)[:n]

    bias_estimate = match_peaks_to_expected(expected_compounds, top_peaks, rt_tolerance=rt_tolerance_bias_estimation)
    signed_deltas = [row['observed_rt'] - row['expected_rt'] for row in bias_estimate if row['observed_rt'] is not None]
    rt_bias = statistics.median(signed_deltas) if signed_deltas else 0.0
    if signed_deltas:
        print(f"  Estimated RT bias {rt_bias:+.2f}s from {len(signed_deltas)}/{n} most intense peaks")

    return match_peaks_to_expected(expected_compounds, observed_peaks, rt_tolerance=rt_tolerance, rt_bias=rt_bias)

def compute_run_metrics(results: list[MatchResult]) -> RunMetrics:
    """
    Summarize one file's matches: the quartiles (Q1, Q2, Q3) of delta RT, S/N and
    FWHM, the mean intensity across the matched compounds, and how many of the
    expected compounds were found.

    The mean is None if no compound was matched; a quartile triplet additionally
    requires at least 2 matched compounds (otherwise also None, since a spread
    computed from a single point isn't meaningful).
    """
    matched = [row for row in results if row['observed_rt'] is not None]
    out: dict[str, Any] = {'n_found': len(matched), 'n_expected': len(results)}
    for key in ('delta_rt', 'sn', 'fwhm'):
        values = [row[key] for row in matched]  # type: ignore[literal-required]
        out[f'{key}_quantiles'] = (tuple(statistics.quantiles(values, n=4, method='inclusive'))
                                    if len(values) >= 2 else None)
    intensities = [row['intensity'] for row in matched]
    out['mean_intensity'] = statistics.fmean(intensities) if intensities else None
    return out  # type: ignore[return-value]

def analyze_chromatogram(
    chromatogram: MSChromatogram,
    expected_compounds: dict[str, ExpectedCompound],
) -> FileAnalysis:
    """
    Pick the peaks of a chromatogram, match them against the expected compounds
    and summarize the matches -- the whole analysis, without any file I/O.

    Args:
        chromatogram: the chromatogram to analyze (typically a TIC)
        expected_compounds: dict mapping compound name -> ExpectedCompound (RT in
            seconds), as returned by load_expected_compounds

    Returns:
        A FileAnalysis. If no peak could be picked at all, its `matches` are empty
        (there is nothing to match against) and `metrics` report 0 of 0 found.
    """
    smoothed_chrom, observed_peaks, noise_levels = pick_chromatographic_peaks(chromatogram)
    matches = match_with_rt_bias_correction(expected_compounds, observed_peaks) if observed_peaks else []
    return FileAnalysis(
        chromatogram=chromatogram,
        smoothed_chromatogram=smoothed_chrom,
        noise_levels=noise_levels,
        observed_peaks=observed_peaks,
        matches=matches,
        metrics=compute_run_metrics(matches),
    )

def analyze_mzml(filename: str | Path, expected_compounds: dict[str, ExpectedCompound]) -> FileAnalysis:
    """
    Read an mzML file, compute its TIC and analyze it (see analyze_chromatogram).

    Raises:
        ValueError: if the file's TIC is empty (nothing to analyze).
    """
    chromatogram = compute_chromatogram(read_mzml_file(filename))
    if chromatogram.size() == 0:
        raise ValueError(f"Empty chromatogram in {filename}")
    return analyze_chromatogram(chromatogram, expected_compounds)

def store_chromatogram_mzml(chromatogram: MSChromatogram, output_filename: str | Path) -> None:
    """Write a single chromatogram to an mzML file."""
    exp = MSExperiment()
    exp.addChromatogram(chromatogram)
    MzMLFile().store(str(output_filename), exp)

def write_results_to_tsv(results: list[MatchResult], output_filename: str | Path) -> None:
    """Write results to a TSV file."""
    with open(output_filename, 'w') as f:
        f.write("compound_name\tRT\tobserved_RT\tdelta_RT\tintensity\tFWHM\tSN\n")
        for row in results:
            f.write(f"{row['compound_name']}\t"
                   f"{row['expected_rt']}\t"
                   f"{row['observed_rt'] if row['observed_rt'] is not None else ''}\t"
                   f"{row['delta_rt'] if row['delta_rt'] is not None else ''}\t"
                   f"{row['intensity'] if row['intensity'] is not None else ''}\t"
                   f"{row['fwhm'] if row['fwhm'] is not None else ''}\t"
                   f"{row['sn'] if row['sn'] is not None else ''}\n")

def write_observed_peaks_to_tsv(observed_peaks: list[PeakInfo], output_filename: str | Path) -> None:
    """Write all observed peaks (rt, intensity, fwhm, S/N) to a TSV file."""
    with open(output_filename, 'w') as f:
        f.write("RT\tintensity\tFWHM\tSN\n")
        for peak in observed_peaks:
            f.write(f"{peak['rt']}\t{peak['intensity']}\t{peak['fwhm']}\t{peak['sn']}\n")

