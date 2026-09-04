#!/usr/bin/env python3
"""
Command line tool to QC GC/MS runs: pick chromatographic peaks in one or more
mzML files, match them against a list of expected (target) compounds and report
how well each run performed.

This module only orchestrates; every step lives in a module of its own that can
also be used on its own:
    gcms_analysis  -- mzML/TIC handling, peak picking, matching, metrics, compound list
    gcms_plotting  -- chromatogram PNGs
    gcms_mzqc      -- mzQC report

Each input mzML file is processed independently and gets its own set of outputs,
named after it and written to OUT_DIR (see -o):
  <name>.mzQC                     run-level and per-compound metrics as mzQC
  <name>.chromatogram_peaks.tsv   one row per expected compound: RT, observed RT,
                                  delta RT, intensity, FWHM, S/N
  <name>.observed_peaks.tsv       all picked peaks: RT, intensity, FWHM, S/N
  <name>.chromatogram_plot.png    chromatogram with picked/matched peaks
  <name>.chromatogram_plot_logY.png   the same with a log-scaled y-axis
  smoothed/<name>.TIC_smoothed.mzML   the smoothed TIC the peaks were picked on

Usage:
    python GC_MS_pickPeaks_and_match.py [-c COMPOUNDS] [-o OUT_DIR] file1.mzML [file2.mzML ...]
    python GC_MS_pickPeaks_and_match.py --write-compounds-template

    -compounds/-c is the CSV list of expected compounds (default: compounds.csv);
    --write-compounds-template writes an adaptable example list to that path.
    -o/--out-dir is interpreted relative to the directory of each input mzML file
    (an absolute path is used as-is) and defaults to './results'.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gcms_analysis import (
    FileAnalysis,
    analyze_mzml,
    load_expected_compounds,
    store_chromatogram_mzml,
    write_compounds_template,
    write_observed_peaks_to_tsv,
    write_results_to_tsv,
)
from gcms_mzqc import write_mzqc
from gcms_plotting import plot_chromatogram

DEFAULT_COMPOUNDS_FILE = 'compounds.csv'
DEFAULT_OUT_DIR = './results'

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Pick chromatographic peaks in mzML files and match them against expected compounds.",
        epilog="Example: python %(prog)s -compounds compounds.csv -o results sample1.mzML sample2.mzML")
    parser.add_argument('mzml_files', nargs='*', metavar='mzML_file',
                        help="one or more input mzML files")
    parser.add_argument('-compounds', '-c', '--compounds', default=DEFAULT_COMPOUNDS_FILE,
                        metavar='COMPOUNDS',
                        help="CSV file listing the expected compounds and their retention times "
                             "(columns 'compound_name' plus 'rt_minutes' or 'rt_seconds') "
                             "(default: %(default)s)")
    parser.add_argument('--write-compounds-template', nargs='?', const=True, default=None,
                        metavar='PATH',
                        help="write an example compound list (two compounds, with comments) that "
                             "you can adapt, then exit; PATH defaults to the -compounds file")
    parser.add_argument('-o', '--out-dir', default=DEFAULT_OUT_DIR, metavar='OUT_DIR',
                        help="directory for all output files, relative to the directory of each "
                             "input mzML file (absolute paths are used as-is) "
                             "(default: %(default)s)")
    return parser.parse_args(argv)

def get_output_dir(input_path: Path, out_dir: str) -> Path:
    """Resolve (and create) the output directory for one input file."""
    output_dir = input_path.parent / out_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir

def write_outputs(analysis: FileAnalysis, input_path: Path, output_dir: Path) -> None:
    """
    Write everything this tool produces for a single mzML file into output_dir,
    each output named after the input file: the mzQC, the two TSVs and the linear
    and log-y plot.
    """
    output_base = output_dir / input_path.name

    mzqc_filename = output_base.with_suffix('.mzQC')
    write_mzqc([(input_path, analysis.metrics, analysis.matches)], mzqc_filename)
    print(f"mzQC written to {mzqc_filename}")

    tsv_filename = output_base.with_suffix('.chromatogram_peaks.tsv')
    write_results_to_tsv(analysis.matches, tsv_filename)
    print(f"Results written to {tsv_filename}")

    observed_peaks_filename = output_base.with_suffix('.observed_peaks.tsv')
    write_observed_peaks_to_tsv(analysis.observed_peaks, observed_peaks_filename)
    print(f"Observed peaks written to {observed_peaks_filename}")

    plot_chromatogram(analysis.chromatogram, analysis.smoothed_chromatogram, analysis.noise_levels,
                      analysis.observed_peaks, analysis.matches,
                      output_base.with_suffix('.chromatogram_plot.png'))
    plot_chromatogram(analysis.chromatogram, analysis.smoothed_chromatogram, analysis.noise_levels,
                      analysis.observed_peaks, analysis.matches,
                      output_base.with_suffix('.chromatogram_plot_logY.png'), log_y=True)

def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.write_compounds_template is not None:
        template_path = Path(args.compounds if args.write_compounds_template is True
                             else args.write_compounds_template)
        write_compounds_template(template_path)
        print(f"Compound list template written to {template_path} -- adapt it and rerun with "
              f"-compounds {template_path}")
        return 0

    if not args.mzml_files:
        print("No input mzML file given. See --help.")
        return 1

    try:
        expected_compounds = load_expected_compounds(args.compounds)
    except FileNotFoundError:
        print(f"Compound list not found: {args.compounds}\n"
              f"Write an adaptable example list with: "
              f"python {Path(__file__).name} --write-compounds-template {args.compounds}")
        return 1
    except ValueError as e:
        print(f"Invalid compound list: {e}")
        return 1
    print(f"Read {len(expected_compounds)} expected compounds from {args.compounds}")

    for input_file in args.mzml_files:
        input_path = Path(input_file)

        if not input_path.exists():
            print(f"File not found: {input_file}")
            continue

        print(f"Processing {input_file}...")
        try:
            # All outputs of this file go here; '<out_dir>/<input name>.<...>'
            output_dir = get_output_dir(input_path, args.out_dir)

            analysis = analyze_mzml(input_path, expected_compounds)

            # The smoothed TIC is written even if nothing was picked on it, since
            # that is exactly the case one wants to look at afterwards
            smoothed_dir = output_dir / "smoothed"
            smoothed_dir.mkdir(parents=True, exist_ok=True)
            smoothed_TIC_filename = (smoothed_dir / input_path.name).with_suffix('.TIC_smoothed.mzML')
            store_chromatogram_mzml(analysis.smoothed_chromatogram, smoothed_TIC_filename)
            print(f"Smoothed chromatogram written to {smoothed_TIC_filename}")

            if not analysis.observed_peaks:
                print(f"No peaks detected in {input_file}. Skipping.")
                continue

            print(f"  Matched {analysis.metrics['n_found']}/{analysis.metrics['n_expected']} expected compounds")
            write_outputs(analysis, input_path, output_dir)

        except Exception as e:
            print(f"Error processing {input_file}: {e}")

    return 0

if __name__ == "__main__":
    sys.exit(main())
