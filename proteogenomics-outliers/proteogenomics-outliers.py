import datetime
import json
import math
import os
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.special
import seaborn as sns
from scipy.spatial.distance import cdist

#settings 

MIN_VAR, MIN_CORR, K, LMBD = 0.01, 0.9, 20, 3
SCALING = "robust"
METRIC = "euclidean"
NEIGHBOUR_MODE = "drop_self"
FLAG_LOW, FLAG_HIGH = 0.33, 0.66
ATTR_YLIM = 1.0
PCA_SEED = 123          

QUAMETER_XLSX = "Mtb-120-outlier-metrics.xlsx"
PTXQC = {
    "noMBR": ("report_v1.1.6_txt_defaultParam_heatmap__proteogenomics-outliers.txt",
              "report_v1.1.6_txt_defaultParam_filename_sort_proteogenomics-outliers.txt"),
    "MBR": ("report_v1.1.6_txt_withMBR_heatmap_proteogenomics-outliers.txt",
            "report_v1.1.6_txt_withMBR_filename_sort_proteogenomics-outliers.txt"),
}
RAW_LOCATION = "ftp://massive.ucsd.edu/v01/MSV000081205/raw/"

#preprocessing + loOp


def preprocess(data, min_variance, min_corr, scaling_mode):
    variances = data.var(axis=0, ddof=0)
    data = data[data.columns[(variances > min_variance).values]]

    mask = np.ones(len(data.columns), dtype=bool)
    corr = data.corr().abs()
    for i in range(len(data.columns)):
        if not mask[i]:
            continue
        for j in range(i + 1, len(data.columns)):
            if corr.iloc[i, j] > min_corr:
                mask[j] = False
    data = data[data.columns[mask]]

    if scaling_mode == "robust":
        centre = data.median()
        scale = data.quantile(0.75) - data.quantile(0.25)
    else:
        centre = data.mean()
        scale = data.std(ddof=0)
    return (data - centre) / scale


def detect_outliers_loop(data, n_neighbors, lmbd=3, metric="euclidean",
                         neighbour_mode="drop_self"):
    X = data.to_numpy(dtype=float)
    first = 1 if neighbour_mode == "drop_self" else 2

    dist = cdist(X, X, metric={"manhattan": "cityblock"}.get(metric, metric))
    order = np.argsort(dist, axis=1, kind="stable")
    knn_ind = order[:, first:first + n_neighbors]
    knn_dists = np.take_along_axis(dist, knn_ind, axis=1)

    pdists = np.sqrt(np.sum(knn_dists ** 2, axis=1) / n_neighbors)
    plofs = pdists / pdists[knn_ind].mean(axis=1) - 1
    nplof = lmbd * np.sqrt(np.mean(plofs ** 2))

    loop = np.clip(scipy.special.erf(plofs / (nplof * np.sqrt(2))), 0, 1)
    loop = pd.Series(loop, index=data.index)

    attributions = pd.DataFrame(index=data.index, columns=data.columns, dtype=float)
    for i in range(X.shape[0]):
        per_feat = ((X[[i]] - X[knn_ind[i]]) ** 2).mean(axis=0)
        attributions.iloc[i] = per_feat / per_feat.sum()
    return loop, attributions


def run_loop(data):
    d = preprocess(data, MIN_VAR, MIN_CORR, SCALING)
    return detect_outliers_loop(d, K, LMBD, METRIC, NEIGHBOUR_MODE)

#quameter tsv  ->  mzQC

QUAMETER_METRICS = [
    ("MS:4000050", "XIC50 fraction", ["XIC-WideFrac"]),
    ("MS:4000051", "XIC-FWHM quantiles", ["XIC-FWHM-Q1", "XIC-FWHM-Q2", "XIC-FWHM-Q3"]),
    ("MS:4000182", "XIC-Height quartile ratios", ["XIC-Height-Q2", "XIC-Height-Q3", "XIC-Height-Q4"]),
    ("MS:4000053", "chromatography duration", ["RT-Duration"]),
    ("MS:4000183", "TIC quarters RT fraction", ["RT-TIC-Q1", "RT-TIC-Q2", "RT-TIC-Q3", "RT-TIC-Q4"]),
    ("MS:4000184", "MS1 quarter RT fraction", ["RT-MS-Q1", "RT-MS-Q2", "RT-MS-Q3", "RT-MS-Q4"]),
    ("MS:4000185", "MS2 quarter RT fraction", ["RT-MSMS-Q1", "RT-MSMS-Q2", "RT-MSMS-Q3", "RT-MSMS-Q4"]),
    ("MS:4000186", "MS1 TIC-change quartile ratios", ["MS1-TIC-Change-Q2", "MS1-TIC-Change-Q3", "MS1-TIC-Change-Q4"]),
    ("MS:4000187", "MS1 TIC quartile ratios", ["MS1-TIC-Q2", "MS1-TIC-Q3", "MS1-TIC-Q4"]),
    ("MS:4000059", "number of MS1 spectra", ["MS1-Count"]),
    ("MS:4000065", "fastest frequency for MS level 1 collection", ["MS1-Freq-Max"]),
    ("MS:4000061", "MS1 density quantiles", ["MS1-Density-Q1", "MS1-Density-Q2", "MS1-Density-Q3"]),
    ("MS:4000060", "number of MS2 spectra", ["MS2-Count"]),
    ("MS:4000066", "fastest frequency for MS level 2 collection", ["MS2-Freq-Max"]),
    ("MS:4000062", "MS2 density quantiles", ["MS2-Density-Q1", "MS2-Density-Q2", "MS2-Density-Q3"]),
]
CHARGE_TABLES = [
    ("MS:4000063", "MS2 known precursor charges fractions",
     ["1", "2", "3", "4", "5", "6"],
     ["MS2-PrecZ-1", "MS2-PrecZ-2", "MS2-PrecZ-3", "MS2-PrecZ-4", "MS2-PrecZ-5", "MS2-PrecZ-more"]),
    ("MS:4000064", "MS2 unknown and likely precursor charges fractions",
     ["1", "2"], ["MS2-PrecZ-likely-1", "MS2-PrecZ-likely-multi"]),
]

RAW_FORMAT = {"accession": "MS:1000563", "name": "Thermo RAW format",
              "description": "Thermo Scientific RAW file format."}
TSV_FORMAT = {"accession": "MS:1000914", "name": "tab delimited text format",
              "description": "A file format that has two or more columns of "
                             "tabular data where each column is separated by a TAB character."}
PSI_MS_CV = {"name": "Proteomics Standards Initiative Mass Spectrometry Ontology",
             "uri": "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/"
                    "1b696a6c2cc6cb0df1639e3b0d0342d41152aa2c/psi-ms.obo",
             "version": "4.1.259"}


def py(v):
    return v.item() if hasattr(v, "item") else v


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def raw_input_file(name, time):
    return {"name": name, "location": f"{RAW_LOCATION}{name}.raw",
            "fileFormat": RAW_FORMAT,
            "fileProperties": [{
                "accession": "MS:1000747", "name": "completion time",
                "description": "The time that a data processing action was finished.",
                "value": time}]}


def write_quameter_mzqc(tsv, path):
    runs = []
    for _, row in tsv.iterrows():
        name = row["Filename"].removesuffix(".raw")
        qms = []
        for acc, mname, cols in QUAMETER_METRICS:
            vals = [py(row[c]) for c in cols]
            qms.append({"accession": acc, "name": mname,
                        "value": vals[0] if len(vals) == 1 else vals})
        for acc, mname, charges, cols in CHARGE_TABLES:
            qms.append({"accession": acc, "name": mname,
                        "value": {"MS:1000041": charges,
                                  "UO:0000191": [py(row[c]) for c in cols]}})
        runs.append({
            "metadata": {
                "label": name,
                "inputFiles": [raw_input_file(name, row["StartTimeStamp"])],
                "analysisSoftware": [{
                    "accession": "MS:1009002", "name": "QuaMeter IDFree",
                    "version": "1.1.21142",
                    "uri": "http://proteowizard.sourceforge.net/"}]},
            "qualityMetrics": qms})
    doc = {"mzQC": {"version": "1.0.0", "creationDate": now(),
                    "description": "QuaMeter IDFree metrics for 120 Mtb runs.",
                    "controlledVocabularies": [PSI_MS_CV],
                    "runQualities": runs}}
    with open(path, "w") as f:
        json.dump(doc, f, indent=2)
    print(f"  wrote {path}")


def read_quameter_mzqc(path):
    mz = json.load(open(path))["mzQC"]
    rows, names, times = {}, {}, {}
    for run in mz["runQualities"]:
        lab = run["metadata"]["label"]
        rows[lab] = {}
        times[lab] = run["metadata"]["inputFiles"][0]["fileProperties"][0]["value"]
        for q in run["qualityMetrics"]:
            v, a = q["value"], q["accession"]
            if isinstance(v, bool):
                continue
            if isinstance(v, (int, float)):
                rows[lab][a] = v
                names[a] = q["name"]
            elif isinstance(v, list):
                for i, x in enumerate(v):
                    rows[lab][f"{a}_{i}"] = x
                    names[f"{a}_{i}"] = f"{q['name']} Q{i + 1}"
    return pd.DataFrame(rows).T, names, times, mz["runQualities"]

#ptxqc heatmaps


def load_ptxqc(input_dir, heatmap, filename_sort):
    heat = pd.read_csv(os.path.join(input_dir, heatmap), sep="\t", quotechar='"',
                       escapechar="\\").set_index("fc.raw.file")
    names = pd.read_csv(os.path.join(input_dir, filename_sort), sep="\t", comment="#")
    mapping = dict(zip(names["new.Name"], names["orig.Name"]))
    missing = set(heat.index) - set(mapping)
    if missing:
        raise ValueError(f"no real name for: {sorted(missing)[:5]}")
    heat.index = [mapping[i] for i in heat.index]
    return heat


def write_loop_mzqc(scores, path, label, description, source_file, run_times,
                    ptxqc_version=None):
    runs = list(run_times)
    scores = scores.reindex(runs)
    software = []
    if ptxqc_version:
        software.append({
            "accession": "MS:1003162", "name": "PTX-QC",
            "description": "Proteomics (PTX) - QualityControl (QC) software for "
                           "QC report generation and visualization.",
            "version": ptxqc_version, "uri": "https://github.com/cbielow/PTXQC"})
    software.append({
        "accession": "MS:1000799", "name": "custom unreleased software tool",
        "description": "A software tool that has not yet been released. The value "
                       "should describe the software. Please do not use this term "
                       "for publicly available software - contact the PSI-MS working "
                       "group in order to have another CV term added.",
        "value": "quality control LoOP outlier detection", "version": "2016-03-14",
        "uri": "https://pubs.acs.org/doi/full/10.1021/acs.jproteome.6b00028"})

    src = os.path.abspath(source_file)
    input_files = [raw_input_file(r, run_times[r]) for r in runs]
    input_files.append({"name": os.path.splitext(os.path.basename(src))[0],
                        "location": "file://" + src, "fileFormat": TSV_FORMAT})

    doc = {"mzQC": {
        "version": "1.0.0", "creationDate": now(), "description": description,
        "controlledVocabularies": [
            PSI_MS_CV,
            {"name": "NCI Thesaurus OBO Edition",
             "uri": "http://purl.obolibrary.org/obo/ncit.owl"}],
        "setQualities": [{
            "metadata": {"label": label, "inputFiles": input_files,
                         "analysisSoftware": software},
            "qualityMetrics": [{
                "accession": "MS:4000201", "name": "run outlier score LoOP",
                "description": "Outlier score for MS runs based on the local outlier "
                               "probabilities (LoOP) algorithm applied to QC metrics "
                               "from all runs.",
                "value": {"MS:4000086": runs,
                          "NCIT:C54154": [round(float(s), 5) for s in scores]}}]}]}}
    with open(path, "w") as f:
        json.dump(doc, f, indent=2)
    print(f"  wrote {path}")

#robust pca 

R_PCA = r"""
suppressMessages(library(MASS))
args <- commandArgs(trailingOnly = TRUE)
MtbDF <- read.csv(args[1], row.names = 1, check.names = FALSE)
set.seed(as.integer(args[3]))
robust.cov <- cov.rob(MtbDF)
robust.cor.1 <- robust.cov
robust.cor.1$cov <- cov2cor(robust.cov$cov)
metrics.pca.1 <- sweep(MtbDF, 2, sqrt(diag(robust.cov$cov)), FUN = "/")
robust.cor.1$center <- robust.cov$center / sqrt(diag(robust.cov$cov))
pc.cr <- princomp(metrics.pca.1, covmat = robust.cor.1, scores = TRUE)
distMatrix <- as.matrix(dist(pc.cr$scores[, 1:5], method = "euclidean"))
Medians <- apply(distMatrix, 1, median, na.rm = TRUE)
write.csv(data.frame(run = rownames(MtbDF), median = Medians), args[2],
          row.names = FALSE)
"""


def pca_features(runs):
    def get(run, acc):
        return next(q["value"] for q in run["qualityMetrics"] if q["accession"] == acc)
    rows = {}
    for run in runs:
        rows[run["metadata"]["label"]] = {
            "XIC50": get(run, "MS:4000050"),
            "XICFWHMQ2": get(run, "MS:4000051")[1],
            "XICHeightQ2": get(run, "MS:4000182")[1],
            "TICQuartRTMiddle": sum(get(run, "MS:4000183")[1:3]),
            "MS1QuartRTMiddle": sum(get(run, "MS:4000184")[1:3]),
            "MS2QuartRTMiddle": sum(get(run, "MS:4000185")[1:3]),
            "MS1TICDeltaQuartRatioQ2": get(run, "MS:4000186")[1],
            "MS1TICQuartRatioQ2": get(run, "MS:4000187")[1],
            "MS1Count": get(run, "MS:4000059"),
            "MS2Count": get(run, "MS:4000060"),
            "MS1DensityQ2": get(run, "MS:4000061")[1],
            "MS2DensityQ2": get(run, "MS:4000062")[1],
            "MS2PreZ2": get(run, "MS:4000063")["UO:0000191"][1],
            "MS2PreZ3": get(run, "MS:4000063")["UO:0000191"][2],
            "MS1FastRate": get(run, "MS:4000065"),
            "MS2FastRate": get(run, "MS:4000066"),
        }
    return pd.DataFrame(rows).T


def robust_pca_score(runs):
    with tempfile.TemporaryDirectory() as tmp:
        feat, med, script = (os.path.join(tmp, f) for f in
                             ("features.csv", "medians.csv", "pca.R"))
        pca_features(runs).to_csv(feat, float_format="%.17g")
        with open(script, "w") as f:
            f.write(R_PCA)
        subprocess.run(["Rscript", script, feat, med, str(PCA_SEED)], check=True)
        m = pd.read_csv(med, index_col=0)["median"]
    return (m - m.min()) / (m.max() - m.min())

#figure

plt.style.use(["seaborn-v0_8-white", "seaborn-v0_8-paper"])
sns.set_context("paper")

PARALLEL_COLS = ["QuaMeter PCA", "QuaMeter LoOP", "PTXQC no MBR"]
DISPLAY_RENAME = {"PTXQC no MBR": "PTXQC LoOP"}


def parallel_plot(df, flagged, fname, width=4.5):
    d = df[PARALLEL_COLS].rename(columns=DISPLAY_RENAME)
    d["flagged"] = flagged
    fig, ax = plt.subplots(figsize=(width, width / 1.618))
    pd.plotting.parallel_coordinates(
        d.sort_values("flagged"), "flagged", ax=ax,
        color=["lightgray", "#6da7de", "#d82222"], axvlines=True, clip_on=False)
    ax.set_ylabel("Outlier score", fontsize=14)
    ax.tick_params(axis="both", labelsize=12.5)
    ax.get_legend().remove()
    plt.tight_layout()
    plt.savefig(fname, dpi=2200, bbox_inches="tight")
    plt.close()
    print(f"  wrote {fname}")


X_TICK, Y_TICK, TITLE, YLABEL = 16.5, 16, 17, 18


def split_label_into_two_lines(label):
    label = str(label)
    words = label.split()
    if len(words) < 2:
        m = max(1, len(label) // 2)
        return f"{label[:m]}\n{label[m:]}"
    at = min(range(1, len(words)),
             key=lambda i: abs(len(" ".join(words[:i])) - len(" ".join(words[i:]))))
    return f"{' '.join(words[:at])}\n{' '.join(words[at:])}"


def attribution_panel(ax, attr_row, score, run, rename, ylim):
    top = attr_row.nlargest(5)
    idx = [rename.get(c, c) for c in top.index]
    ax.bar(range(5), top.values, color="#6da7de", width=0.88)
    ax.set_xticks(range(5))
    ax.set_xticklabels([split_label_into_two_lines(l) for l in idx],
                       rotation=90, ha="center", va="top")
    ax.tick_params(axis="x", pad=3, labelsize=X_TICK)
    ax.tick_params(axis="y", labelsize=Y_TICK, length=5, width=1.3)
    ax.margins(x=0.02)
    ax.set_ylim(0, ylim)
    ax.set_title(f"{run} (LoOP score = {score:.2f})", fontsize=TITLE)
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(1.3)
        sp.set_color("black")


def attribution_figure(panels, fname, ncols=2, ylim=0.5, panel_w=5.0):
    nrows = math.ceil(len(panels) / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharey=True,
                             figsize=(panel_w * ncols, 4.8 * nrows))
    axes = np.atleast_2d(axes)
    for (row, score, run, rename), ax in zip(panels, axes.flatten()):
        attribution_panel(ax, row, score, run, rename, ylim)
    for ax in axes.flatten()[len(panels):]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        if ax.get_visible():
            ax.set_ylabel("Feature attribution ratio", fontsize=YLABEL)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.28, wspace=0.08, hspace=0.55)
    plt.savefig(fname, dpi=600, bbox_inches="tight", transparent=True)
    plt.close()
    print(f"  wrote {fname}")

#main


def main():
    input_dir = sys.argv[1] if len(sys.argv) > 1 else "../data/input"
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "../data/output"
    os.makedirs(output_dir, exist_ok=True)
    out = lambda f: os.path.join(output_dir, f)

    #quaMeter metrics -> mzQC
    print("QuaMeter metrics -> mzQC")
    tsv = pd.read_excel(os.path.join(input_dir, QUAMETER_XLSX))
    qm_mzqc = out("Mtb-120-outlier-metrics.mzQC")
    write_quameter_mzqc(tsv, qm_mzqc)
    metrics, names, run_times, qm_runs = read_quameter_mzqc(qm_mzqc)

    #loop scores
    print("LoOP outlier scores")
    qm_score, qm_attr = run_loop(metrics)
    write_loop_mzqc(qm_score, out("Mtb-120-outlier-loop-quameter.mzQC"),
                    "Mtb-120-outlier-loop-quameter",
                    "LoOP outlier scores for 120 Mtb runs from QuaMeter IDFree metrics.",
                    qm_mzqc, run_times)

    ptx_score = {}
    for cond, mbr in (("noMBR", "disabled"), ("MBR", "enabled")):
        heatmap, filename_sort = PTXQC[cond]
        heat = load_ptxqc(input_dir, heatmap, filename_sort)
        ptx_score[cond], _ = run_loop(heat)
        write_loop_mzqc(ptx_score[cond], out(f"Mtb-120-outlier-loop-ptxqc-{cond}.mzQC"),
                        f"Mtb-120-outlier-loop-ptxqc-{cond}",
                        "LoOP outlier scores for 120 Mtb runs from PTX-QC v1.1.6 "
                        f"metrics (MaxQuant, match-between-runs {mbr}).",
                        os.path.join(input_dir, heatmap), run_times,
                        ptxqc_version="1.1.6")

    print("Robust PCA (R)")
    pca = robust_pca_score(qm_runs)

    print("Figures")
    df = pd.DataFrame({"QuaMeter LoOP": qm_score, "QuaMeter PCA": pca,
                       "PTXQC no MBR": ptx_score["noMBR"]}).reindex(list(run_times))
    m = df[PARALLEL_COLS].max(axis=1)
    flagged = (m > FLAG_HIGH).astype(int) + (m > FLAG_LOW)
    print("  flag counts:", flagged.value_counts().sort_index().to_dict())
    parallel_plot(df, flagged, out("outlier_scores.png"))

    s = qm_score.sort_values(ascending=False)
    top = s.index[s > FLAG_HIGH]
    print("  QuaMeter LoOP > %.2f: " % FLAG_HIGH
          + ", ".join(f"{i} = {s[i]:.3f}" for i in top))
    attribution_figure([(qm_attr.loc[i], s[i], i, names) for i in top],
                       out("attributions_quameter.png"), ylim=ATTR_YLIM)


if __name__ == "__main__":
    main()
