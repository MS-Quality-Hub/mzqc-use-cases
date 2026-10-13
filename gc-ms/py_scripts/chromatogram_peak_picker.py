#!/usr/bin/env python3
"""
Custom chromatographic peak picker, adapted from OpenMS's PeakPickerChromatogram
(src/openms/source/ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp).

PeakPickerChromatogram computes a per-point signal-to-noise (S/N) estimate
internally (via SignalToNoiseEstimatorMedian) purely to decide how far to
extend each peak's left/right border, then discards it. This module
reimplements the same "corrected" picking algorithm on top of pyOpenMS's
existing building blocks (GaussFilter/SavitzkyGolayFilter, PeakPickerHiRes,
SignalToNoiseEstimatorMedianChrom) so that the apex S/N can be captured and
reported for every picked peak.
"""

from __future__ import annotations

from typing import Any, TypedDict
from pyopenms import (
    MSChromatogram,
    GaussFilter,
    SavitzkyGolayFilter,
    PeakPickerHiRes,
    SignalToNoiseEstimatorMedianChrom,
)

class PickedPeak(TypedDict):
    """One picked chromatographic peak with its apex S/N and integration bounds."""
    rt: float
    intensity: float
    fwhm: float | None
    sn: float
    left_rt: float
    right_rt: float
    integrated_intensity: float

DEFAULT_PARAMS: dict[str, Any] = {
    'use_gauss': True,
    'gauss_width': 2.0,             # Gaussian width in seconds, estimated peak size
    'sgolay_frame_length': 15,
    'sgolay_polynomial_order': 3,
    'signal_to_noise': 10.0,        # Minimum S/N ratio for peak detection/extension
    'sn_win_len': 100.0,           # S/N estimator window length (seconds)
    'sn_bin_count': 300,             # S/N estimator histogram bin count
    'peak_width': -1.0,             # Force a minimal peak width (seconds); -1 disables
}


def _smooth(chromatogram: MSChromatogram, cfg: dict[str, Any]) -> MSChromatogram:
    """Smooth a copy of the chromatogram, mirroring PeakPickerChromatogram's smoothing step."""
    smoothed = MSChromatogram(chromatogram)
    if cfg['use_gauss']:
        gauss = GaussFilter()
        params = gauss.getParameters()
        params.setValue("gaussian_width", cfg['gauss_width'])
        gauss.setParameters(params)
        gauss.filter(smoothed)
    else:
        sgolay = SavitzkyGolayFilter()
        params = sgolay.getParameters()
        params.setValue("frame_length", cfg['sgolay_frame_length'])
        params.setValue("polynomial_order", cfg['sgolay_polynomial_order'])
        sgolay.setParameters(params)
        sgolay.filter(smoothed)
    return smoothed


def _pick_seeds(smoothed_chrom: MSChromatogram, cfg: dict[str, Any]) -> MSChromatogram:
    """
    Find seed peaks (apex RT, intensity, FWHM) via PeakPickerHiRes, as PeakPickerChromatogram does.

    Note: PeakPickerChromatogram's internal PeakPickerHiRes instance has its signal_to_noise
    fixed at the class default of 1.0 at construction time and is never updated when the
    user later sets "signal_to_noise" on the outer object (see PeakPickerChromatogram.cpp,
    constructor vs. updateMembers_()) -- that setting only affects the boundary-extension
    stage below. Reproduced here verbatim rather than "fixed", to stay faithful to the
    reference algorithm's actual (if surprising) behavior.
    """
    picker = PeakPickerHiRes()
    params = picker.getParameters()
    params.setValue("signal_to_noise", 1.0)
    params.setValue("report_FWHM", "true")
    params.setValue("report_FWHM_unit", "absolute")
    picker.setParameters(params)

    picked = MSChromatogram()
    # check_spacings=False: not applicable to chromatograms (equal RT spacing assumed)
    picker.pick(smoothed_chrom, picked, False)
    return picked


def _find_closest_index(chromatogram: MSChromatogram, target_rt: float, start_idx: int = 0) -> int:
    """Mirrors PeakPickerChromatogram::findClosestPeak_: closest index to target_rt, searching from start_idx onward."""
    n = chromatogram.size()
    i = start_idx
    while i < n:
        if target_rt < chromatogram[i].getRT():
            if i > 0 and abs(target_rt - chromatogram[i - 1].getRT()) < abs(target_rt - chromatogram[i].getRT()):
                i -= 1
            return i
        i += 1
    return i


def pick_chromatographic_peaks(
    chromatogram: MSChromatogram, **overrides: Any
) -> list[MSChromatogram | list[PickedPeak] | list[float]]:
    """
    Pick chromatographic peaks, replicating PeakPickerChromatogram's "corrected" algorithm,
    additionally recording the apex signal-to-noise ratio for every picked peak.

    Returns:
        Triple of a smoothed chromatogram, a list of picked peaks as dictionaries with keys
        'rt', 'intensity', 'fwhm', 'sn', 'left_rt', 'right_rt', 'integrated_intensity', and the
        per-point noise level (same length/RTs as the smoothed chromatogram; intensity / S/N,
        with NaN where S/N is 0, i.e. too few points fell into that window).
    """
    cfg: dict[str, Any] = {**DEFAULT_PARAMS, **overrides}

    smoothed_chrom = _smooth(chromatogram, cfg)
    picked_peaks = _pick_seeds(smoothed_chrom, cfg)

    fwhm_arrays = picked_peaks.getFloatDataArrays()
    fwhm_data = fwhm_arrays[0] if len(fwhm_arrays) > 0 else None
    if fwhm_data is not None and len(fwhm_data) != picked_peaks.size():
        raise ValueError("FWHM data does not match picked peaks")

    n_points = smoothed_chrom.size()
    snt = SignalToNoiseEstimatorMedianChrom()
    snt_params = snt.getParameters()
    snt_params.setValue("win_len", cfg['sn_win_len'])
    snt_params.setValue("bin_count", cfg['sn_bin_count'])
    snt.setParameters(snt_params)
    snt.init(smoothed_chrom)

    noise_levels: list[float] = []
    for i in range(n_points):
        sn = snt.getSignalToNoise(i)
        noise_levels.append(smoothed_chrom[i].getIntensity() / sn if sn > 0.0 else float('nan'))

    signal_to_noise = cfg['signal_to_noise']
    peak_width = cfg['peak_width']

    print(f"Found {picked_peaks.size()} picked peaks")

    results: list[PickedPeak] = []
    current_idx = 0
    for i in range(picked_peaks.size()):
        central_rt = picked_peaks[i].getRT()
        current_idx = _find_closest_index(smoothed_chrom, central_rt, current_idx)
        apex_idx = current_idx
        apex_sn = snt.getSignalToNoise(apex_idx)

        # extend left, mirroring PeakPickerChromatogram::pickChromatogram_
        k = 2
        while ((apex_idx - k + 1) > 0
               and (smoothed_chrom[apex_idx - k].getIntensity() < smoothed_chrom[apex_idx - k + 1].getIntensity()
                    or (peak_width > 0.0 and abs(smoothed_chrom[apex_idx - k].getRT() - central_rt) < peak_width))
               and (signal_to_noise <= 0.0 or snt.getSignalToNoise(apex_idx - k) >= signal_to_noise)):
            k += 1
        left_idx = apex_idx - k + 1

        # extend right
        k = 2
        while ((apex_idx + k) < n_points
               and (smoothed_chrom[apex_idx + k].getIntensity() < smoothed_chrom[apex_idx + k - 1].getIntensity()
                    or (peak_width > 0.0 and abs(smoothed_chrom[apex_idx + k].getRT() - central_rt) < peak_width))
               and (signal_to_noise <= 0.0 or snt.getSignalToNoise(apex_idx + k) >= signal_to_noise)):
            k += 1
        right_idx = apex_idx + k - 1

        # for peak integration, use the raw (unsmoothed) chromatogram
        integrated_intensity = sum(chromatogram[j].getIntensity() for j in range(left_idx, right_idx + 1))

        results.append({
            'rt': central_rt,
            'intensity': picked_peaks[i].getIntensity(),
            'fwhm': fwhm_data[i] if fwhm_data is not None else None,
            'sn': apex_sn,
            'left_rt': chromatogram[left_idx].getRT(),
            'right_rt': chromatogram[right_idx].getRT(),
            'integrated_intensity': integrated_intensity,
        })

    return [smoothed_chrom, results, noise_levels]
