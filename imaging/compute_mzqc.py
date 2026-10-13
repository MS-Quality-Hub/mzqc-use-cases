"""Compute mass spectrometry imaging QC metrics and write them as an mzQC file.

    MS:4000228  total ion current image                  matrix
    MS:4000229  total ion current image quantiles        n-tuple (Q1, Q2, Q3)
    MS:4000230  base peak fraction image                 matrix
    MS:4000231  base peak fraction image quantiles       n-tuple (Q1, Q2, Q3)
    MS:4000232  reference m/z error image                matrix
    MS:4000233  reference m/z error image quantiles      n-tuple (Q1, Q2, Q3)
    MS:4000234  acquired pixel fraction                  single value
    MS:4000236  duplicate pixel coordinate count         single value
    MS:4000237  pixel grid dimensions                    n-tuple (X, Y)
    MS:4000235  QC reference m/z                         metadata.cvParameters

Two existing psi-ms metrics are emitted alongside them:

    MS:4000099  number of empty MS1 scans
    MS:4000069  m/z acquisition range

Usage
-----
    python compute_mzqc.py DATA.imzML --ref-mz 1296.6848 --window-ppm 100 \
        --spectrum-type profile -o ../data/output/DATA.mzQC
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pyopenms as oms
from mzqc.MZQCFile import (
    AnalysisSoftware,
    ControlledVocabulary,
    CvParameter,
    InputFile,
    JsonSerialisable,
    MetaDataParameters,
    MzQcFile,
    QualityMetric,
    RunQuality,
)

PSI_MS_VERSION = "4.2.4"
PSI_MS_URI = (
    "https://github.com/HUPO-PSI/psi-ms-CV/releases/download/"
    f"v{PSI_MS_VERSION}/psi-ms.obo"
)

# accession -> (name, description, unit accession | None)
TERMS: dict[str, tuple[str, str, str | None]] = {
    "MS:1003611": ("imzML format", "IMSIS imzML imaging file format.", None),
    "MS:1003430": (
        "OpenMS", "A software for LC-MS data management and analysis.", None,
    ),
    "MS:4000235": (
        "QC reference m/z",
        "The theoretical m/z of the mass reference ion used by a mass "
        "reference based QC metric.",
        "MS:1000040",
    ),
    "MS:1001412": (
        "search tolerance plus value",
        "Search tolerance above the theoretical value, given in relative or "
        "absolute units.",
        "UO:0000169",
    ),
    "MS:1001413": (
        "search tolerance minus value",
        "Search tolerance below the theoretical value, given in relative or "
        "absolute units.",
        "UO:0000169",
    ),
    "MS:1000035": (
        "peak picking",
        "Spectral peak processing conducted on the acquired data to convert "
        "profile data to centroided data.",
        None,
    ),
    "MS:4000228": (
        "total ion current image",
        "The total ion current (MS:1000285) of each spectrum, arranged as a "
        "matrix in which the row and column indices correspond to the y and "
        "x pixel coordinates of the acquisition grid. Positions that were "
        "not acquired are reported as null.",
        "MS:1000043",
    ),
    "MS:4000229": (
        "total ion current image quantiles",
        "The first to n-th quantile of the total ion current (MS:1000285) "
        "values over all acquired spectra of an imaging acquisition. A "
        "value triplet represents the quartiles Q1, Q2, Q3. The number of "
        "values in the tuple implies the quantile mode.",
        "MS:1000043",
    ),
    "MS:4000230": (
        "base peak fraction image",
        "The base peak intensity (MS:1000505) of each centroid spectrum "
        "(MS:1000127) divided by the total ion current (MS:1000285) of the "
        "same spectrum, arranged as a matrix in which the row and column "
        "indices correspond to the y and x pixel coordinates of the "
        "acquisition grid. Positions that were not acquired are reported as "
        "null.",
        "UO:0000191",
    ),
    "MS:4000231": (
        "base peak fraction image quantiles",
        "The first to n-th quantile of the base peak fraction image "
        "(MS:4000230) values over all acquired centroid spectra of an "
        "imaging acquisition. A value triplet represents the quartiles Q1, "
        "Q2, Q3. The number of values in the tuple implies the quantile "
        "mode.",
        "UO:0000191",
    ),
    "MS:4000232": (
        "reference m/z error image",
        "The relative deviation in parts per million between the observed "
        "m/z of a mass reference ion and its theoretical m/z, calculated as "
        "1E6 x (observed m/z - theoretical m/z) / theoretical m/z, for each "
        "spectrum, arranged as a matrix in which the row and column indices "
        "correspond to the y and x pixel coordinates of the acquisition "
        "grid. The theoretical m/z of the reference MUST be reported using "
        "QC reference m/z (MS:4000235). Positions that were not acquired, "
        "and positions in which the reference ion was not detected, are "
        "reported as null.",
        "UO:0000169",
    ),
    "MS:4000233": (
        "reference m/z error image quantiles",
        "The first to n-th quantile of the reference m/z error values over "
        "all spectra of an imaging acquisition in which the mass reference "
        "ion was detected. A value triplet represents the quartiles Q1, Q2, "
        "Q3. The number of values in the tuple implies the quantile mode.",
        "UO:0000169",
    ),
    "MS:4000234": (
        "acquired pixel fraction",
        "The number of spectra with unique pixel coordinates in an imaging "
        "acquisition divided by the number of positions of the smallest "
        "rectangular grid that encloses all acquired pixel coordinates.",
        "UO:0000191",
    ),
    "MS:4000236": (
        "duplicate pixel coordinate count",
        "The number of distinct x/y pixel coordinates in a singular "
        "imaging acquisition at which more than one spectrum of the same "
        "type was recorded.",
        "UO:0000189",
    ),
    "MS:4000237": (
        "pixel grid dimensions",
        "The dimensions, in pixels, of the smallest rectangular grid that "
        "encloses all acquired pixel coordinates of an imaging acquisition. "
        "A value pair (X, Y) represents the extent along x and y, in that "
        "order.",
        "UO:0000189",
    ),
    "MS:4000099": (
        "number of empty MS1 scans",
        "Number of MS1 scans where the scans' peaks intensity sums to 0 "
        "(i.e. no peaks or only 0-intensity peaks).",
        "UO:0000189",
    ),
    "MS:4000069": (
        "m/z acquisition range",
        "Upper and lower limit of m/z precursor values at which MSn spectra "
        "are recorded.",
        "MS:1000040",
    ),
}

UNIT_NAMES: dict[str, str] = {
    "MS:1000040": "m/z",
    "MS:1000043": "intensity unit",
    "UO:0000191": "fraction",
    "UO:0000169": "parts per million",
    "UO:0000189": "count unit",
}


# ==========================================================================
# metric computation
# ==========================================================================
def pick_peaks(mie: oms.MSImagingExperiment) -> oms.MSImagingExperiment:
    """Centroid every spectrum with PeakPickerHiRes (S/N = 0), keeping the geometry."""
    picker = oms.PeakPickerHiRes()
    params = picker.getParameters()
    params.setValue("signal_to_noise", 0.0)
    picker.setParameters(params)

    picked = oms.MSImagingExperiment(picker.pickExperiment(mie.getMSExperiment(), True))
    picked.setGeometry(mie.getGeometry())
    return picked


def duplicate_pixel_coordinates(mie: oms.MSImagingExperiment) -> int:
    """Count of distinct x/y pixel coordinates recorded by more than one spectrum of one type.

    Reads each spectrum's own ``imzml:x``/``imzml:y`` MetaValue rather than the
    imaging geometry, which keeps only the first spectrum per coordinate."""
    def key(s: oms.MSSpectrum) -> tuple:
        return (
            int(s.getMetaValue("imzml:x")),
            int(s.getMetaValue("imzml:y")),
            s.getMSLevel(),
            s.getInstrumentSettings().getPolarity(),
            tuple(p.getMZ() for p in s.getPrecursors()),
        )

    return sum(1 for n in Counter(map(key, mie.getMSExperiment().getSpectra())).values() if n > 1)


def quantiles(image: np.ndarray) -> list[float]:
    """Q1, Q2, Q3 over the finite entries of a matrix metric."""
    return [float(q) for q in np.percentile(image[np.isfinite(image)], [25, 50, 75])]


def compute(
    imzml_path: Path, ref_mz: float, window_ppm: float, spectrum_type: str,
) -> dict[str, Any]:
    """Load an imzML and compute every metric in a single pass over the pixels.

    ``spectrum_type`` is ``"profile"`` (peak-picked here first) or ``"centroid"``
    (used as it stands). Returns the metric values keyed by accession.
    """
    mie = oms.MSImagingExperiment()
    oms.ImzMLFile().load(str(imzml_path), mie)
    n_duplicates = duplicate_pixel_coordinates(mie)

    if spectrum_type == "profile":
        mie = pick_peaks(mie)
    exp = mie.getMSExperiment()

    # Grid = bounding box of the pixels in the geometry (one spectrum per coordinate).
    pixels = mie.getGeometry().pixels_struct()
    x0, y0 = int(pixels["x"].min()), int(pixels["y"].min())
    width = int(pixels["x"].max()) - x0 + 1
    height = int(pixels["y"].max()) - y0 + 1

    tic = np.full((height, width), np.nan)
    bpf = np.full((height, width), np.nan)
    mzerr = np.full((height, width), np.nan)

    spectra = exp.getSpectra()
    tol = ref_mz * window_ppm / 1e6
    n_empty = 0
    for px in pixels:
        yy, xx = int(px["y"]) - y0, int(px["x"]) - x0
        s = spectra[int(px["spectrum_index"])]
        total = s.calculateTIC()
        tic[yy, xx] = total
        if total > 0:
            bpf[yy, xx] = s.get_base_peak()[1] / total
        else:
            n_empty += 1
        # search for reference m/z as highest peak inside window defined by ppm tolerance
        k = s.findHighestInWindow(ref_mz, tol, tol)
        if k >= 0:
            mzerr[yy, xx] = (s[k].getMZ() - ref_mz) / ref_mz * 1e6

    exp.updateRanges()
    return {
        "MS:4000228": tic,
        "MS:4000229": quantiles(tic),
        "MS:4000230": bpf,
        "MS:4000231": quantiles(bpf),
        "MS:4000232": mzerr,
        "MS:4000233": quantiles(mzerr),
        "MS:4000234": len(pixels) / (width * height),
        "MS:4000236": n_duplicates,
        "MS:4000237": [width, height],
        "MS:4000099": n_empty,
        "MS:4000069": [float(exp.getMinMZ()), float(exp.getMaxMZ())],
    }


# ==========================================================================
# mzQC assembly
# ==========================================================================
def _unit(accession: str | None) -> CvParameter | None:
    return None if accession is None else CvParameter(accession=accession, name=UNIT_NAMES[accession])


def _cv(accession: str, value: Any = None) -> CvParameter:
    name, description, unit = TERMS[accession]
    return CvParameter(
        accession=accession, name=name, description=description, value=value, unit=_unit(unit),
    )


def _metric(accession: str, value: Any) -> QualityMetric:
    """A QualityMetric; matrices go out as nested lists with null where the value is absent."""
    name, description, unit = TERMS[accession]
    if isinstance(value, np.ndarray):
        value = [[float(v) if np.isfinite(v) else None for v in row] for row in value]
    return QualityMetric(
        accession=accession, name=name, description=description, value=value, unit=_unit(unit),
    )


def to_mzqc(
    imzml_path: Path, metrics: dict[str, Any], ref_mz: float, window_ppm: float,
    spectrum_type: str,
) -> MzQcFile:
    """One runQuality for one imzML/ibd pair."""
    # imzML is a two-file format; both halves carry MS:1003611.
    input_files = [
        InputFile(name=p.name, location=p.resolve().as_uri(), fileFormat=_cv("MS:1003611"))
        for p in (imzml_path, imzml_path.with_suffix(".ibd"))
    ]
    openms_name, openms_def, _ = TERMS["MS:1003430"]
    software = [
        AnalysisSoftware(
            accession="MS:1003430",
            name=openms_name,
            description=openms_def,
            value=f"pyOpenMS {oms.__version__}",
            version=oms.__version__,
            uri="https://openms.de",
        )
    ]
    cv_parameters = [
        _cv("MS:4000235", value=ref_mz),
        _cv("MS:1001412", value=window_ppm),
        _cv("MS:1001413", value=window_ppm),
    ]
    if spectrum_type == "profile":
        cv_parameters.append(
            _cv(
                "MS:1000035",
                value="PeakPickerHiRes, signal_to_noise=0; all per-pixel metrics "
                "computed on the centroided spectra",
            )
        )

    run = RunQuality(
        metadata=MetaDataParameters(
            label=imzml_path.stem,
            inputFiles=input_files,
            analysisSoftware=software,
            cvParameters=cv_parameters,
        ),
        qualityMetrics=[_metric(acc, value) for acc, value in metrics.items()],
    )
    return MzQcFile(
        version="1.0.0",
        creationDate=datetime.now(),
        description=(
            "Imaging QC metrics for the mzQC usecase manuscript."
        ),
        runQualities=[run],
        controlledVocabularies=[
            ControlledVocabulary(
                name="Proteomics Standards Initiative Mass Spectrometry Ontology",
                uri=PSI_MS_URI,
                version=PSI_MS_VERSION,
            ),
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("imzml", type=Path, help="imzML file (its .ibd must sit beside it)")
    parser.add_argument("--ref-mz", type=float, required=True,
                        help="theoretical m/z of the mass reference ion (MS:4000235)")
    parser.add_argument("--window-ppm", type=float, required=True,
                        help="reference search window, +/- ppm (no default, see above)")
    parser.add_argument("--spectrum-type", choices=("profile", "centroid"), required=True,
                        help="representation of the input spectra; profile is peak-picked first")
    parser.add_argument("-o", "--out", type=Path, default=None,
                        help="output mzQC path (default: out/<imzML stem>.mzQC)")
    args = parser.parse_args()

    metrics = compute(args.imzml, args.ref_mz, args.window_ppm, args.spectrum_type)
    mzqc = to_mzqc(args.imzml, metrics, args.ref_mz, args.window_ppm, args.spectrum_type)
    out = args.out or Path("out") / f"{args.imzml.stem}.mzQC"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(JsonSerialisable.to_json(mzqc, readability=1))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
