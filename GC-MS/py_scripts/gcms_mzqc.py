#!/usr/bin/env python3
"""
mzQC export of the per-file GC/MS QC metrics (via pymzqc).

Self-contained and usable on its own with the metrics of gcms_analysis:

    from gcms_analysis import analyze_mzml, load_expected_compounds
    from gcms_mzqc import write_mzqc

    analysis = analyze_mzml('sample.mzML', load_expected_compounds('compounds.csv'))
    write_mzqc([('sample.mzML', analysis.metrics, analysis.matches)],
               'result_stats_summary.mzQC')
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, NamedTuple, TYPE_CHECKING

if TYPE_CHECKING:
    from gcms_analysis import MatchResult, RunMetrics

# CV terms for the per-file metrics, in output order. MS:4000216-4000219 are custom
# terms (see psi-ms.obo, appended at the end of the official CV, pending submission
# to PSI-MS) since the PSI-MS CV has no existing "across targeted/ID-based compounds"
# metric for retention time deviation, S/N, FWHM, or intensity -- only ID-free,
# whole-spectrum variants (e.g. MS:4000122 "MS1 signal-to-noise ratio mean").
#
# The three distributions are reported as quartile triplets (Q1, Q2, Q3), following
# the n-tuple pattern of e.g. MS:4000051 "XIC-FWHM quantiles", rather than as
# separate median and IQR metrics: Q2 is the median and Q3 - Q1 the IQR.
class CvTerm(NamedTuple):
    """A CV term used as a QualityMetric, with the CV term of its unit (if any)."""
    accession: str
    name: str
    unit: tuple[str, str] | None = None  # (accession, name) of the unit term

_SECOND = ('UO:0000010', 'second')
_INTENSITY = ('MS:1000043', 'intensity unit')
TARGETED_COMPOUND_CV_TERMS: dict[str, CvTerm] = {
    'delta_rt_quantiles': CvTerm('MS:4000216', 'targeted compound retention time deviation quantiles', _SECOND),
    'sn_quantiles':       CvTerm('MS:4000217', 'targeted compound signal-to-noise ratio quantiles'),
    'fwhm_quantiles':     CvTerm('MS:4000218', 'targeted compound peak width at half-maximum quantiles', _SECOND),
    'mean_intensity':     CvTerm('MS:4000219', 'targeted compound intensity mean', _INTENSITY),
}

# CV terms for the per-compound table metrics: the metric term itself plus the CV
# terms naming its columns (mzQC identifies table columns by accession, not by name).
# The compound names are data, not vocabulary, so any compound list can be reported.
COMPOUND_NAME_COLUMN = ('MS:4000222', 'targeted compound name')
TARGETED_COMPOUND_TABLE_CV_TERMS: dict[str, tuple[CvTerm, CvTerm]] = {
    # match key -> (table metric term, value column term)
    'delta_rt':  (CvTerm('MS:4000220', 'targeted compound retention time deviations'),
                  CvTerm('MS:4000223', 'targeted compound retention time deviation', _SECOND)),
    'intensity': (CvTerm('MS:4000221', 'targeted compound intensities'),
                  CvTerm('MS:4000224', 'targeted compound intensity', _INTENSITY)),
}

class SoftwareInfo(NamedTuple):
    """The software that computed the metrics, for the mzQC 'analysisSoftware' metadata."""
    name: str            # written as the CV term's value, i.e. what the software is called
    version: str         # required by the mzQC schema
    uri: str = ''        # optional; omitted from the file when empty
    accession: str = 'MS:1000799'
    cv_name: str = 'custom unreleased software tool'

DEFAULT_SOFTWARE = SoftwareInfo(name='GC_MS_pickPeaks_and_match.py', version='1.0.0')

def _table_metrics(matches: list[MatchResult], mzqc: Any) -> list[Any]:
    """
    Build the per-compound table QualityMetrics (delta RT and intensity) from one
    file's matches. Compounds without a matched peak are left out of the table
    entirely -- per the CV term definitions, a missing row means 'not detected'.
    """
    quality_metrics: list[Any] = []
    for key, (table_term, value_term) in TARGETED_COMPOUND_TABLE_CV_TERMS.items():
        rows = [row for row in matches if row[key] is not None]  # type: ignore[literal-required]
        if not rows:
            continue
        quality_metrics.append(mzqc.QualityMetric(
            accession=table_term.accession,
            name=table_term.name,
            value={COMPOUND_NAME_COLUMN[0]: [row['compound_name'] for row in rows],
                   value_term.accession: [row[key] for row in rows]},  # type: ignore[literal-required]
        ))
    return quality_metrics

def write_mzqc(file_metrics: list[tuple[str | Path, RunMetrics, list[MatchResult]]],
                output_filename: str | Path,
                software: SoftwareInfo = DEFAULT_SOFTWARE) -> None:
    """
    Write an mzQC file (via pymzqc) with one RunQuality per processed mzML file,
    reporting per file:
      - the run-level aggregates of TARGETED_COMPOUND_CV_TERMS: the quartiles of
        delta RT, S/N and FWHM (as n-tuple metrics) and the mean intensity across
        the matched compounds. Aggregates that are None (no compound matched in
        that file, or fewer than 2 for the quartiles) are omitted rather than
        written as null.
      - the per-compound tables of TARGETED_COMPOUND_TABLE_CV_TERMS: delta RT and
        intensity per compound, with the compound names in their own column.

    Args:
        file_metrics: list of (input_path, metrics, matches) tuples, where `metrics`
            is the dict returned by compute_run_metrics (keys
            'delta_rt_quantiles', 'sn_quantiles', 'fwhm_quantiles', 'mean_intensity')
            and `matches` are that file's per-compound MatchResults
        software: the tool that computed the metrics, reported as the required
            'analysisSoftware' metadata of every run
    """
    import json
    from mzqc import MZQCFile as mzqc

    cv = mzqc.ControlledVocabulary(
        name='PSI-MS',
        uri='https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo',
        version='4.1.257',
    )

    analysis_software = mzqc.AnalysisSoftware(
        accession=software.accession, name=software.cv_name, value=software.name,
        version=software.version, uri=software.uri,
    )

    run_qualities: list[Any] = []
    for input_path, metrics, matches in file_metrics:
        input_file = mzqc.InputFile(
            location=Path(input_path).resolve().as_uri(),
            name=Path(input_path).name,
            fileFormat=mzqc.CvParameter(accession='MS:1000584', name='mzML format'),
        )
        # pymzqc materializes every optional list as [], which the mzQC schema rejects
        # ('should be non-empty'); setting them to None makes the serializer drop them
        input_file.fileProperties = None
        quality_metrics: list[Any] = []
        for key, term in TARGETED_COMPOUND_CV_TERMS.items():
            value = metrics[key]  # type: ignore[literal-required]
            if value is None:
                continue
            # n-tuple metrics go out as a JSON array, single values as a scalar
            if isinstance(value, tuple):
                value = list(value)
            kwargs: dict[str, Any] = dict(accession=term.accession, name=term.name, value=value)
            if term.unit:
                # the unit has to be a CV parameter, not a bare accession, or the
                # file fails semantic validation
                kwargs['unit'] = mzqc.CvParameter(accession=term.unit[0], name=term.unit[1])
            quality_metrics.append(mzqc.QualityMetric(**kwargs))
        quality_metrics.extend(_table_metrics(matches, mzqc))
        metadata = mzqc.MetaDataParameters(
            label=Path(input_path).name,  # required, and unique across the runs
            inputFiles=[input_file],
            analysisSoftware=[analysis_software],
        )
        metadata.cvParameters = None
        run_qualities.append(mzqc.RunQuality(metadata=metadata, qualityMetrics=quality_metrics))

    mzqc_file = mzqc.MzQcFile(
        contactName='GC_MS_pickPeaks_and_match.py',
        description=(
            'Quartiles (Q1, Q2, Q3) of delta RT, S/N and FWHM, the mean intensity, and '
            'the per-compound delta RT and intensity tables, across targeted/expected '
            'GC-QC compounds, per input mzML file. Uses custom CV terms MS:4000216-'
            'MS:4000224, appended to the copy of psi-ms.obo in this project\'s analysis/ '
            'folder pending official PSI-MS submission.'
        ),
        runQualities=run_qualities,
        controlledVocabularies=[cv],
    )
    mzqc_file.setQualities = None  # empty list again -- this file has no set qualities

    with open(output_filename, 'w') as f:
        json.dump({'mzQC': mzqc_file}, f, default=mzqc.JsonSerialisable.complex_handler,
                   cls=mzqc.MzqcJSONEncoder, indent=2)
