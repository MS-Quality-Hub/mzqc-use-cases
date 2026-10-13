from __future__ import annotations
import argparse
import datetime
import gzip
import json
import math
import os
import re
import sys
import urllib.parse
from collections import Counter
from dataclasses import dataclass, field
from typing import Optional
import pandas as pd

DEFAULT_INPUT = "./metabolomics_repository_qc.csv"
DEFAULT_OBO = "./psi-ms.obo"
DEFAULT_OUTDIR = "./metabolomics_mzQC"

DEFAULT_ACCESSION_BASE = 4000216

DISTINCT_MZ_DECIMALS = 2

FLOAT_DECIMALS = 4

SOFTWARE_NAME = "metabolomics repository QC pipeline"
SOFTWARE_VERSION = "1.0.0"
SOFTWARE_URI = "https://github.com/HUPO-PSI/mzQC"

CV_NAME = "Proteomics Standards Initiative Mass Spectrometry Ontology"

PLACEHOLDER_ACCESSION = "MS:4000XXX"

DESCRIBE_PROPOSED = False
ACCESSION_RE = re.compile(r"^[A-Z]+:[A-Z0-9]+$")



@dataclass
class ProposedTerm:
    key: str                       
    name: str
    definition: str
    value_parent: str              
    value_type: str                
    unit_accession: str
    unit_name: str
    categories: list              
    columns: list                 
    comment: Optional[str] = None
    value_concept: Optional[tuple] = None   
    accession: str = PLACEHOLDER_ACCESSION  


CAT_ID_FREE = ("MS:4000009", "ID free metric")
CAT_SINGLE_RUN = ("MS:4000012", "single run based metric")
CAT_MS = ("MS:4000019", "MS metric")
CAT_MS1 = ("MS:4000021", "MS1 metric")
CAT_MS2 = ("MS:4000022", "MS2 metric")
CAT_ISOLATION = ("MS:1000792", "isolation window attribute")

UNIT_COUNT = ("UO:0000189", "count unit")
UNIT_MZ = ("MS:1000040", "m/z")

VALUE_PARENT_SINGLE = "MS:4000003"
VALUE_PARENT_NTUPLE = "MS:4000004"

_D = DISTINCT_MZ_DECIMALS

PROPOSED_TERMS = [
    ProposedTerm(
        key="ms1_pos",
        name="number of positive polarity MS1 spectra",
        definition="The number of MS1 spectra acquired in positive polarity "
                   "in a single MS run.",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:int",
        unit_accession=UNIT_COUNT[0], unit_name=UNIT_COUNT[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS1],
        columns=["MS1_pos_count"],
    ),
    ProposedTerm(
        key="ms1_neg",
        name="number of negative polarity MS1 spectra",
        definition="The number of MS1 spectra acquired in negative polarity "
                   "in a single MS run.",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:int",
        unit_accession=UNIT_COUNT[0], unit_name=UNIT_COUNT[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS1],
        columns=["MS1_neg_count"],
    ),
    ProposedTerm(
        key="ms2_pos",
        name="number of positive polarity MS2 spectra",
        definition="The number of MS2 spectra acquired in positive polarity "
                   "in a single MS run.",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:int",
        unit_accession=UNIT_COUNT[0], unit_name=UNIT_COUNT[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS2],
        columns=["MS2_pos_count"],
    ),
    ProposedTerm(
        key="ms2_neg",
        name="number of negative polarity MS2 spectra",
        definition="The number of MS2 spectra acquired in negative polarity "
                   "in a single MS run.",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:int",
        unit_accession=UNIT_COUNT[0], unit_name=UNIT_COUNT[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS2],
        columns=["MS2_neg_count"],
    ),
    ProposedTerm(
        key="ms3_pos",
        name="number of positive polarity MS3+ spectra",
        definition="The number of MS3 and higher (MSn, n>2) spectra acquired "
                   "in positive polarity in a single MS run.",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:int",
        unit_accession=UNIT_COUNT[0], unit_name=UNIT_COUNT[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS],
        columns=["MS3+_pos_count"],
    ),
    ProposedTerm(
        key="ms3_neg",
        name="number of negative polarity MS3+ spectra",
        definition="The number of MS3 and higher (MSn, n>2) spectra acquired "
                   "in negative polarity in a single MS run.",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:int",
        unit_accession=UNIT_COUNT[0], unit_name=UNIT_COUNT[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS],
        columns=["MS3+_neg_count"],
    ),
    ProposedTerm(
        key="ms1_observed_bounds",
        name="MS1 observed m/z bounds",
        definition="The minimum and maximum m/z among detected peaks in MS1 "
                   "spectra in a single MS run, reported as a [min, max] tuple.",
        value_parent=VALUE_PARENT_NTUPLE, value_type="xsd:float",
        unit_accession=UNIT_MZ[0], unit_name=UNIT_MZ[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS1],
        columns=["MS1_Min_MZ", "MS1_Max_MZ"],
        value_concept=("STATO:0000035", "range"),
    ),
    ProposedTerm(
        key="ms2_range_above_prec",
        name="MS2 median acquisition range above precursor",
        definition="The median m/z distance between the MS2 precursor m/z and "
                   "the upper limit of its fragment scan window, taken across "
                   "all MS2 spectra in a single MS run.",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:float",
        unit_accession=UNIT_MZ[0], unit_name=UNIT_MZ[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS2],
        columns=["MS2_Aquisition_Range_Above_Prec"],
    ),
    ProposedTerm(
        key="ms2_isolation_width",
        name="MS2 precursor isolation window median width",
        definition="The median width of the precursor isolation window applied "
                   "during MS2 fragmentation in a single MS run, measured as "
                   "the distance from the lower to the upper isolation window "
                   "bound, that is the sum of the isolation window lower offset "
                   "(MS:1000828) and the isolation window upper offset "
                   "(MS:1000829). The width is independent of the isolation "
                   "window target m/z and remains defined if the target falls "
                   "outside the window bounds.",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:float",
        unit_accession=UNIT_MZ[0], unit_name=UNIT_MZ[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS2, CAT_ISOLATION],
        columns=["MS2_Isolation_Width"],
    ),
    ProposedTerm(
        key="distinct_precursor_mz",
        name="number of distinct precursor m/z values",
        definition=f"The number of distinct selected ion m/z values "
                   f"(MS:1000744) recorded in MS2 scan headers in a single MS "
                   f"run, where two values are considered identical if they "
                   f"agree after rounding to {_D} decimal places.",
        comment="Reads the selected ion m/z (data-dependent selection). "
                "Differs from number of distinct isolation target m/z values, "
                "which reads the isolation window target m/z (MS:1000827).",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:int",
        unit_accession=UNIT_COUNT[0], unit_name=UNIT_COUNT[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS2],
        columns=["Num_Unique_Precursor_MZ"],
    ),
    ProposedTerm(
        key="distinct_target_mz",
        name="number of distinct isolation target m/z values",
        definition=f"The number of distinct isolation window target m/z values "
                   f"(MS:1000827) programmed for fragmentation in a single MS "
                   f"run, where two values are considered identical if they "
                   f"agree after rounding to {_D} decimal places.",
        comment="Reads the isolation window target m/z (the isolation scheme; "
                "fixed windows in DIA/PRM/SRM). See number of distinct "
                "precursor m/z values for the selected ion m/z counterpart; "
                "the two differ when isolation is offset from the selected ion "
                "or when a target is selected repeatedly.",
        value_parent=VALUE_PARENT_SINGLE, value_type="xsd:int",
        unit_accession=UNIT_COUNT[0], unit_name=UNIT_COUNT[1],
        categories=[CAT_ID_FREE, CAT_SINGLE_RUN, CAT_MS2, CAT_ISOLATION],
        columns=["Num_Unique_Target_MZ"],
    ),
]

#metrics that depend on the upstream distinctness rule
DISTINCTNESS_KEYS = {"distinct_precursor_mz", "distinct_target_mz"}

#existing cv terms reused (accession -> fallback name)
EXISTING_TERMS = {
    "MS:4000059": "number of MS1 spectra",
    "MS:4000060": "number of MS2 spectra",
    "MS:4000067": "MS run duration",
}

DERIVED_INPUT_COLS = ["RT_Range_in_Min"]

UNUSED_COLS = [
    "Contains_Profile", "Contains_Centroid",
    "MS1_Min_Acquisition_MZ", "MS1_Max_Acquisition_MZ",
    "MS1_Acquisition_Range", "MS1_MZ_Range",
]


PLACEHOLDER_SUFFIXES = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def assign_accessions(base: int, style: str = "numbered") -> None:
    for i, term in enumerate(PROPOSED_TERMS):
        if style == "shared":
            term.accession = PLACEHOLDER_ACCESSION
        elif style == "unique":
            term.accession = f"MS:4000XX{PLACEHOLDER_SUFFIXES[i]}"
        else:
            term.accession = f"MS:{base + i}"


@dataclass
class CvTerm:
    accession: str
    name: str = ""
    definition: str = ""
    is_a: list = field(default_factory=list)
    units: Optional[str] = None
    value_type: Optional[str] = None


def parse_obo(path_or_text: str, is_text: bool = False):
    if is_text:
        text = path_or_text
    else:
        with open(path_or_text, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()

    m = re.search(r"^data-version:\s*(.+)$", text, re.M)
    version = m.group(1).strip() if m else "unknown"

    terms = {}
    for block in text.split("\n[")[0:]:
        if not block.startswith("Term]"):
            continue
        body = block.split("]", 1)[1]
        acc = name = defn = None
        is_a, units, vtype = [], None, None
        for line in body.splitlines():
            line = line.rstrip()
            if line.startswith("id: "):
                acc = line[4:].strip()
            elif line.startswith("name: "):
                name = line[6:].strip()
            elif line.startswith("def: "):
                dm = re.match(r'def:\s*"(.*)"\s*\[', line)
                defn = dm.group(1) if dm else line[5:].strip()
            elif line.startswith("is_a: "):
                is_a.append(line[6:].split("!")[0].strip())
            elif line.startswith("relationship: has_units "):
                units = line[len("relationship: has_units "):].split("!")[0].strip()
            elif line.startswith("relationship: has_value_type "):
                vtype = line[len("relationship: has_value_type "):].split("!")[0].strip()
        if acc:
            terms[acc] = CvTerm(acc, name or "", defn or "", is_a, units, vtype)
    return version, terms


def load_cv(obo_path: Optional[str]):
    cv_url = "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"
    if obo_path and os.path.exists(obo_path):
        print(f"Loading CV from local file: {obo_path}")
        version, terms = parse_obo(obo_path)
    else:
        print(f"Fetching CV from: {cv_url}")
        import urllib.request
        with urllib.request.urlopen(cv_url) as resp:
            version, terms = parse_obo(resp.read().decode("utf-8", "replace"),
                                       is_text=True)
    print(f"  Loaded {len(terms):,} CV terms (data-version: {version})")
    return version, terms


def check_terms_against_cv(cv_terms, cv_version) -> list:
    """verify every proposed term matches the released cv and adopt the cv wording"""
    problems, notes = [], []
    for term in PROPOSED_TERMS:
        cv = cv_terms.get(term.accession)
        if cv is None:
            problems.append(f"{term.accession} ({term.name}) is not in the CV "
                            f"(data-version {cv_version})")
            continue
        if cv.name != term.name:
            notes.append(f"{term.accession}: name '{term.name}' -> "
                         f"'{cv.name}' (from CV)")
            term.name = cv.name
        if cv.definition and cv.definition != term.definition:
            notes.append(f"{term.accession}: definition updated from CV")
            term.definition = cv.definition
        if cv.value_type and cv.value_type != term.value_type:
            problems.append(f"{term.accession}: value type {term.value_type} "
                            f"but CV says {cv.value_type}")
        if cv.units and cv.units != term.unit_accession:
            problems.append(f"{term.accession}: unit {term.unit_accession} "
                            f"but CV says {cv.units}")
        if cv.is_a and term.value_parent not in cv.is_a:
            problems.append(f"{term.accession}: value parent "
                            f"{term.value_parent} but CV says {cv.is_a}")
    if notes:
        print("  CV is authoritative, adopted its wording:")
        for n in notes:
            print(f"    {n}")
    if problems:
        print("\n  term mismatch refusing to write files:")
        for p_ in problems:
            print(f"    {p_}")
    return problems


def build_cv_uri(cv_version: str) -> str:
    return ("https://github.com/HUPO-PSI/psi-ms-CV/releases/download/"
            f"v{cv_version}/psi-ms.obo")



def render_obo_draft() -> str:
    lines = [
        "! Proposed QC CV terms for HUPO-PSI/psi-ms-CV",
        f"! Provisional accessions; distinctness rule: {DISTINCT_MZ_DECIMALS} decimal places",
        "",
    ]
    parent_names = {VALUE_PARENT_SINGLE: "single value",
                    VALUE_PARENT_NTUPLE: "n-tuple"}
    for t in PROPOSED_TERMS:
        lines.append("[Term]")
        lines.append(f"id: {t.accession}")
        lines.append(f"name: {t.name}")
        lines.append(f'def: "{t.definition}" [PSI:MS]')
        if t.comment:
            lines.append(f"comment: {t.comment}")
        lines.append(f"is_a: {t.value_parent} ! {parent_names[t.value_parent]}")
        for cat_acc, cat_name in t.categories:
            lines.append(f"relationship: has_metric_category {cat_acc} ! {cat_name}")
        lines.append(f"relationship: has_value_type {t.value_type} "
                     f"! The allowed value-type for this CV term")
        lines.append(f"relationship: has_units {t.unit_accession} ! {t.unit_name}")
        if t.value_concept:
            lines.append(f"relationship: has_value_concept {t.value_concept[0]} "
                         f"! {t.value_concept[1]}")
        lines.append("")
    return "\n".join(lines)


def as_int(val):
    if val is None or (isinstance(val, float) and math.isnan(val)) or pd.isna(val):
        return None
    try:
        return int(round(float(val)))
    except (TypeError, ValueError):
        return None


def as_float(val, decimals: int = FLOAT_DECIMALS):
    if val is None or pd.isna(val):
        return None
    try:
        v = float(val)
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return round(v, decimals)


def parse_mri(mri: str):
    parts = mri.split(":", 2)
    if len(parts) >= 3:
        dataset, rel_path = parts[1], parts[2]
    else:
        dataset, rel_path = "unknown", mri
    basename = rel_path.rsplit("/", 1)[-1]
    return dataset, rel_path, basename


def strip_extension(name: str) -> str:
    return name.rsplit(".", 1)[0] if "." in name else name


def build_location(dataset: str, rel_path: str, style: str) -> str:
    quoted = urllib.parse.quote(rel_path, safe="/")
    if style == "usi":
        return f"mzspec:{dataset}:{rel_path}"
    if style == "massive" and dataset.startswith("MSV"):
        return f"ftp://massive.ucsd.edu/{dataset}/{quoted}"
    return f"file:///{urllib.parse.quote(dataset, safe='')}/{quoted}"


def detect_file_format(rel_path: str, cv_terms):
    low = rel_path.lower()
    if low.endswith(".mzml"):
        acc = "MS:1000584"
    elif low.endswith(".mzxml"):
        acc = "MS:1000566"
    elif low.endswith(".raw"):
        acc = "MS:1000563"
    elif low.endswith(".mgf"):
        acc = "MS:1001062"
    else:
        acc = "MS:1000560"          
    ct = cv_terms.get(acc)
    return {"accession": acc, "name": ct.name if ct else "mass spectrometer file format"}


def unique_suffix(rel_path: str, used: set) -> str:
    parts = rel_path.split("/")
    for depth in range(1, len(parts) + 1):
        candidate = "/".join(parts[-depth:])
        if candidate not in used:
            used.add(candidate)
            return candidate
    i = 2
    while f"{rel_path}_{i}" in used:
        i += 1
    used.add(f"{rel_path}_{i}")
    return f"{rel_path}_{i}"

#metric construction 


def cv_description(accession: str, cv_terms) -> Optional[str]:
    term = cv_terms.get(accession)
    if term is None or not term.definition:
        return None
    return term.definition


def _existing_metric(accession, value, unit, cv_terms):
    """emit a metric for a term that already exists in the cv name and description are taken from the cv bc pymzqc errors on any divergence from the ontology"""
    term = cv_terms.get(accession)
    metric = {"accession": accession,
              "name": term.name if term else EXISTING_TERMS[accession]}
    desc = cv_description(accession, cv_terms)
    if desc:
        metric["description"] = desc
    metric["value"] = value
    metric["unit"] = {"accession": unit[0], "name": unit[1]}
    return metric


def build_quality_metrics(row, cv_terms, stats: Counter,
                          drop_implausible: bool = False,
                          strict_distinctness: bool = False):
    metrics = []

    for term in PROPOSED_TERMS:
        if strict_distinctness and term.key in DISTINCTNESS_KEYS:
            continue

        if term.value_parent == VALUE_PARENT_NTUPLE:
            lo = as_float(row.get(term.columns[0]))
            hi = as_float(row.get(term.columns[1]))
            if lo is None or hi is None:
                continue
            if lo > hi:
                stats["ntuple_min_gt_max"] += 1
                if drop_implausible:
                    continue
            value = [lo, hi]
        else:
            raw = row.get(term.columns[0])
            value = as_int(raw) if term.value_type == "xsd:int" else as_float(raw)
            if value is None:
                continue
            if value < 0:
                stats[f"negative:{term.key}"] += 1
                if drop_implausible:
                    continue

        metric = {"accession": term.accession, "name": term.name}
        desc = cv_description(term.accession, cv_terms)
        if desc is None and DESCRIBE_PROPOSED:
            desc = term.definition
        if desc:
            metric["description"] = desc
        metric["value"] = value
        metric["unit"] = {"accession": term.unit_accession,
                          "name": term.unit_name}
        metrics.append(metric)
        stats[f"emitted:{term.accession}"] += 1

    #existing cv terms
    ms1 = (as_int(row.get("MS1_pos_count")) or 0) + (as_int(row.get("MS1_neg_count")) or 0)
    ms2 = (as_int(row.get("MS2_pos_count")) or 0) + (as_int(row.get("MS2_neg_count")) or 0)

    if ms1 > 0:
        metrics.append(_existing_metric(
            "MS:4000059", ms1, ("UO:0000189", "count unit"), cv_terms))
        stats["emitted:MS:4000059"] += 1

    if ms2 > 0:
        metrics.append(_existing_metric(
            "MS:4000060", ms2, ("UO:0000189", "count unit"), cv_terms))
        stats["emitted:MS:4000060"] += 1

    #MS:4000067-> ms run duration -> seconds
    rt_min = as_float(row.get("RT_Range_in_Min"), decimals=6)
    if rt_min is not None and rt_min > 0:
        metrics.append(_existing_metric(
            "MS:4000067", round(rt_min * 60.0, 2),
            ("UO:0000010", "second"), cv_terms))
        stats["emitted:MS:4000067"] += 1
    elif rt_min == 0:
        stats["rt_range_zero"] += 1

    if ms2 == 0 and any(m["accession"] == PROPOSED_TERMS[9].accession and m["value"] > 0
                        for m in metrics):
        stats["precursors_without_ms2"] += 1

    return metrics


def build_run_quality(row, cv_terms, used_labels, stats,
                      location_style, software, drop_implausible,
                      strict_distinctness):
    mri = str(row.get("mri", "")).strip()
    dataset, rel_path, basename = parse_mri(mri)

    metrics = build_quality_metrics(row, cv_terms, stats,
                                    drop_implausible=drop_implausible,
                                    strict_distinctness=strict_distinctness)
    if not metrics:
        return None, dataset

    location = build_location(dataset, rel_path, location_style)

    encoded_path = urllib.parse.quote(rel_path, safe="/")
    seen = used_labels.setdefault(dataset, {"names": set(), "labels": set()})
    name = unique_suffix(encoded_path, seen["names"])
    label = strip_extension(urllib.parse.unquote(name))
    if label in seen["labels"]:
        label = urllib.parse.unquote(name)          
        if label in seen["labels"]:
            i = 2
            while f"{label}_{i}" in seen["labels"]:
                i += 1
            label = f"{label}_{i}"
        stats["label_collisions_resolved"] += 1
    seen["labels"].add(label)

    software_param = {
        "accession": "MS:1000799",
        "name": "custom unreleased software tool",
    }
    soft_desc = cv_description("MS:1000799", cv_terms)
    if soft_desc:
        software_param["description"] = soft_desc
    software_param.update({
        "value": software["name"],
        "version": software["version"],
        "uri": software["uri"],
    })

    return {
        "metadata": {
            "label": label,
            "inputFiles": [{
                "name": name,
                "location": location,
                "fileFormat": detect_file_format(rel_path, cv_terms),
            }],
            "analysisSoftware": [software_param],
        },
        "qualityMetrics": metrics,
    }, dataset

#validation

RFC3339_RE = re.compile(
    r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?(Z|[+-]\d{2}:\d{2})$")


def validate_mzqc(doc, cv_terms, proposed_by_acc):
    """Return a list of error strings; empty means valid."""
    errors = []
    root = doc.get("mzQC")
    if root is None:
        return ["missing root element 'mzQC'"]

    if not re.match(r"^\d+\.\d+\.\d+$", str(root.get("version", ""))):
        errors.append("version must be major.minor.patch")
    if not RFC3339_RE.match(str(root.get("creationDate", ""))):
        errors.append("creationDate is not an RFC3339 date-time")

    run_q = root.get("runQualities", [])
    set_q = root.get("setQualities", [])
    if not run_q and not set_q:
        errors.append("at least one of runQualities / setQualities is required")

    cvs = root.get("controlledVocabularies", [])
    if not cvs:
        errors.append("controlledVocabularies is required and must be non-empty")
    for cv in cvs:
        for key in ("name", "uri"):
            if not cv.get(key):
                errors.append(f"controlledVocabulary missing '{key}'")

    known = set(cv_terms) | set(proposed_by_acc)
    labels = set()
    name_to_location = {}

    def check_param(param, where, expect_metric):
        acc = param.get("accession", "")
        if not ACCESSION_RE.match(acc):
            errors.append(f"{where}: accession '{acc}' violates ^[A-Z]+:[A-Z0-9]+$")
            return
        if not param.get("name"):
            errors.append(f"{where}: missing name for {acc}")
        if acc not in known:
            errors.append(f"{where}: accession {acc} not found in the declared CVs")
            return
        refs = proposed_by_acc.get(acc, [])
        if not refs:
            cv_name = cv_terms[acc].name
            if param.get("name") != cv_name:
                errors.append(f"{where}: name '{param.get('name')}' != CV name '{cv_name}'")
            cv_def = cv_terms[acc].definition
            desc = param.get("description")
            if desc is not None and cv_def and desc != cv_def:
                errors.append(f"{where}: description differs from the CV definition for {acc}")
        if not expect_metric:
            return

        #value type + unit checks
        if refs:
            ref = next((r for r in refs
                        if r.name == param.get("name")), refs[0])
            parent, vtype, unit_acc = ref.value_parent, ref.value_type, ref.unit_accession
        else:
            ref = None
            term = cv_terms[acc]
            parent = (VALUE_PARENT_NTUPLE if VALUE_PARENT_NTUPLE in term.is_a
                      else VALUE_PARENT_SINGLE if VALUE_PARENT_SINGLE in term.is_a
                      else None)
            vtype, unit_acc = term.value_type, term.units

        value = param.get("value")
        if parent == VALUE_PARENT_NTUPLE and not isinstance(value, list):
            errors.append(f"{where}: {acc} is an n-tuple but value is not an array")
        if parent == VALUE_PARENT_SINGLE and isinstance(value, (list, dict)):
            errors.append(f"{where}: {acc} is a single value but value is {type(value).__name__}")
        if vtype == "xsd:int":
            flat = value if isinstance(value, list) else [value]
            if any(not isinstance(v, int) or isinstance(v, bool) for v in flat):
                errors.append(f"{where}: {acc} declares xsd:int but value is not integral")
        if unit_acc:
            unit = param.get("unit")
            if not unit:
                errors.append(f"{where}: {acc} requires a unit ({unit_acc}) but none is given")
            elif isinstance(unit, dict) and unit.get("accession") != unit_acc:
                errors.append(f"{where}: {acc} unit {unit.get('accession')} != CV unit {unit_acc}")

    for kind, container in (("runQuality", run_q), ("setQuality", set_q)):
        for i, quality in enumerate(container):
            where = f"{kind}[{i}]"
            meta = quality.get("metadata", {})
            label = meta.get("label")
            if not label:
                errors.append(f"{where}: metadata.label is required")
            elif label in labels:
                errors.append(f"{where}: duplicate label '{label}'")
            else:
                labels.add(label)

            input_files = meta.get("inputFiles", [])
            if not input_files:
                errors.append(f"{where}: at least one inputFile is required")
            locations = set()
            for f_i, in_file in enumerate(input_files):
                fw = f"{where}.inputFiles[{f_i}]"
                if not in_file.get("name"):
                    errors.append(f"{fw}: name is required")
                loc = in_file.get("location", "")
                if not urllib.parse.urlparse(loc).scheme:
                    errors.append(f"{fw}: location '{loc}' is not a URI")
                if loc in locations:
                    errors.append(f"{fw}: duplicate location within {where}")
                locations.add(loc)
                fname = in_file.get("name", "")
                if fname and fname not in loc:
                    errors.append(f"{fw}: name '{fname}' is not contained in "
                                  f"its location (pymzqc consistency check)")
                if name_to_location.setdefault(fname, loc) != loc:
                    errors.append(f"{fw}: inputFile name '{fname}' maps to more "
                                  f"than one location")
                fmt = in_file.get("fileFormat")
                if not fmt:
                    errors.append(f"{fw}: fileFormat is required")
                else:
                    check_param(fmt, fw + ".fileFormat", expect_metric=False)

            software = meta.get("analysisSoftware", [])
            if not software:
                errors.append(f"{where}: at least one analysisSoftware is required")
            for s_i, soft in enumerate(software):
                sw = f"{where}.analysisSoftware[{s_i}]"
                if not soft.get("version"):
                    errors.append(f"{sw}: version is required")
                check_param(soft, sw, expect_metric=False)

            metrics = quality.get("qualityMetrics", [])
            if not metrics:
                errors.append(f"{where}: at least one qualityMetric is required")
            seen = set()
            for m_i, metric in enumerate(metrics):
                mw = f"{where}.qualityMetrics[{m_i}]"
                acc = metric.get("accession")
                if acc in seen:
                    errors.append(f"{mw}: duplicate metric {acc} within {where}")
                seen.add(acc)
                check_param(metric, mw, expect_metric=True)

    return errors


def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def build_proposed_lookup():
    lookup = {}
    for t in PROPOSED_TERMS:
        lookup.setdefault(t.accession, []).append(t)
    return lookup


def validate_directory(directory, cv_terms, limit=None):
    proposed_by_acc = build_proposed_lookup()
    files = sorted(f for f in os.listdir(directory)
                   if f.endswith(".mzqc") or f.endswith(".mzqc.gz"))
    if limit:
        files = files[:limit]
    bad = 0
    print(f"\nValidating {len(files):,} file(s) in {directory}")
    for name in files:
        with open_maybe_gzip(os.path.join(directory, name)) as fh:
            doc = json.load(fh)
        errs = validate_mzqc(doc, cv_terms, proposed_by_acc)
        if errs:
            bad += 1
            print(f"  [FAIL] {name}")
            for e in errs[:8]:
                print(f"      - {e}")
            if len(errs) > 8:
                print(f"      ... {len(errs) - 8} more")
    print(f"  {len(files) - bad:,} valid, {bad:,} invalid")
    return bad == 0


def convert(csv_path, obo_path, output_dir, n_rows, location_style,
            software, gzip_output, drop_implausible, strict_distinctness,
            cv_uri_override=None):
    print("=" * 64)
    print("metabolomics_repository_qc.csv -> mzQC 1.0.0 (one file per dataset)")
    print("=" * 64)

    os.makedirs(output_dir, exist_ok=True)
    cv_version, cv_terms = load_cv(obo_path)
    cv_uri = cv_uri_override or build_cv_uri(cv_version)
    print(f"  CV URI: {cv_uri}")

    problems = check_terms_against_cv(cv_terms, cv_version)
    if problems:
        raise SystemExit(
            "\naborting: the metric definitions in this script disagree with "
            "the CV\nupdate the script, or point --obo_file at a CV release "
            "that contains the terms\n(they were merged in 4.1.261).")

    metric_cols = sorted({c for t in PROPOSED_TERMS for c in t.columns})
    usecols = ["mri"] + metric_cols + DERIVED_INPUT_COLS

    read_kwargs = {"usecols": usecols, "low_memory": False}
    if n_rows > 0:
        read_kwargs["nrows"] = n_rows
    df = pd.read_csv(csv_path, **read_kwargs)
    for col in metric_cols + DERIVED_INPUT_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    before = len(df)
    df = df.drop_duplicates(subset="mri", keep="first")
    dropped_dupes = before - len(df)
    print(f"\nLoaded {before:,} rows ({dropped_dupes:,} duplicate mri dropped)")

    creation_date = (datetime.datetime.now(datetime.timezone.utc)
                     .isoformat(timespec="seconds").replace("+00:00", "Z"))
    print(f"Creation date: {creation_date}")

    stats = Counter()
    skipped = 0

    #group by dataset and one filer per dataset
    df["_dataset"] = df["mri"].str.split(":", n=2).str[1].fillna("unknown")
    df = df.sort_values("_dataset", kind="stable")

    written = total_runs = 0
    for dataset, group in df.groupby("_dataset", sort=False):
        used_labels = {}
        run_qualities = []
        for row in group.drop(columns="_dataset").to_dict("records"):
            run_quality, _ = build_run_quality(
                row, cv_terms, used_labels, stats, location_style, software,
                drop_implausible, strict_distinctness)
            if run_quality is None:
                skipped += 1
                continue
            run_qualities.append(run_quality)
        if not run_qualities:
            continue

        doc = {"mzQC": {
            "version": "1.0.0",
            "creationDate": creation_date,
            "description": (f"ID-free quality control metrics for "
                            f"{len(run_qualities)} run(s) from dataset {dataset}."),
            "controlledVocabularies": [
                {"name": CV_NAME, "uri": cv_uri, "version": cv_version},
            ],
            "runQualities": run_qualities,
        }}
        fname = re.sub(r"[^\w\-.]", "_", dataset)[:200] or "unknown"
        path = os.path.join(output_dir,
                            f"{fname}.mzqc.gz" if gzip_output else f"{fname}.mzqc")
        if gzip_output:
            with gzip.open(path, "wt", encoding="utf-8") as fh:
                json.dump(doc, fh, indent=2)
        else:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(doc, fh, indent=2)
        written += 1
        total_runs += len(run_qualities)
        del doc, run_qualities
        if written % 500 == 0:
            print(f"  ... {written:,} datasets written")

    obo_path_out = os.path.join(output_dir, "proposed_qc_terms.obo")
    with open(obo_path_out, "w", encoding="utf-8") as fh:
        fh.write(render_obo_draft())

    print(f"\n{'=' * 64}")
    print("DONE")
    print(f"  Datasets (files):    {written:,}")
    print(f"  runQualities:        {total_runs:,}")
    print(f"  Rows with no metric: {skipped:,}")
    print(f"  Output directory:    {output_dir}")
    print(f"  OBO draft:           {obo_path_out}")

    print("\nMetric emission counts")
    for term in PROPOSED_TERMS:
        print(f"  {term.accession}  {term.name:<52} "
              f"{stats[f'emitted:{term.accession}']:>9,}")
    for acc, name in EXISTING_TERMS.items():
        print(f"  {acc}  {name + ' (existing)':<52} {stats[f'emitted:{acc}']:>9,}")


    if not strict_distinctness:
        print("\nNOTE: num_unique_precursor_mz / num_unique_target_mz were computed "
              f"upstream;\n      the {DISTINCT_MZ_DECIMALS} decimal distinctness rule in the term "
              "definitions cannot be\n      applied retroactively to a count. "
              "rerun the extraction at "
              f"{DISTINCT_MZ_DECIMALS} decimals\n      or use --strict_distinctness to omit both metrics.")
    return cv_terms


def parse_arguments():
    p = argparse.ArgumentParser(
        description="Convert metabolomics_repository_qc.csv to mzQC 1.0.0 files.")
    p.add_argument("--input_file", default=DEFAULT_INPUT)
    p.add_argument("--obo_file", default=DEFAULT_OBO)
    p.add_argument("--output_dir", default=DEFAULT_OUTDIR)
    p.add_argument("--n_rows", type=int, default=50,
                   help="rows to process (0 = all)")
    p.add_argument("--accession_base", type=int, default=DEFAULT_ACCESSION_BASE,
                   help="first provisional accession for the proposed terms")
    p.add_argument("--placeholder_style",
                   choices=["shared", "unique", "numbered"], default="numbered",
                   help="numbered = the released accessions from "
                        "--accession_base (default, correct since the terms "
                        "were merged in psi-ms.obo 4.1.261); shared/unique are "
                        "legacy placeholder modes and produce invalid mzQC")
    p.add_argument("--location_style", choices=["file", "massive", "usi"],
                   default="file")
    p.add_argument("--software_name", default=SOFTWARE_NAME)
    p.add_argument("--software_version", default=SOFTWARE_VERSION)
    p.add_argument("--software_uri", default=SOFTWARE_URI)
    p.add_argument("--cv_uri", default=None, help="override the CV URI")
    p.add_argument("--gzip", action="store_true", help="write .mzqc.gz")
    p.add_argument("--drop_implausible", action="store_true",
                   help="skip negative / inverted values instead of emitting them")
    p.add_argument("--strict_distinctness", action="store_true",
                   help="omit the two distinct-m/z metrics")
    p.add_argument("--no_proposed_descriptions", action="store_true",
                   help="omit descriptions for terms not yet in the CV")
    p.add_argument("--validate", action="store_true",
                   help="validate the written files")
    p.add_argument("--validate_limit", type=int, default=0,
                   help="validate only the first N files (0 = all)")
    p.add_argument("--terms_only", action="store_true",
                   help="print the OBO draft and exit")
    return p.parse_args()


def main():
    args = parse_arguments()
    global DESCRIBE_PROPOSED
    DESCRIBE_PROPOSED = not args.no_proposed_descriptions
    assign_accessions(args.accession_base, args.placeholder_style)
    if args.placeholder_style != "numbered":
        print("WARNING: the 11 terms are now released in psi-ms.obo (4.1.261) "
              "as\n         MS:4000216-MS:4000226. Placeholder accessions "
              "produce mzQC that\n         fails semantic validation. Use "
              "--placeholder_style numbered.\n")

    if args.terms_only:
        print(render_obo_draft())
        return

    cv_terms = convert(
        args.input_file, args.obo_file, args.output_dir, args.n_rows,
        args.location_style,
        {"name": args.software_name, "version": args.software_version,
         "uri": args.software_uri},
        args.gzip, args.drop_implausible, args.strict_distinctness,
        cv_uri_override=args.cv_uri)

    if args.validate:
        ok = validate_directory(args.output_dir, cv_terms,
                                limit=args.validate_limit or None)
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
