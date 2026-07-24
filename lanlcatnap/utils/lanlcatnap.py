"""Cross-referencing candidates against known bnAb epitopes (LANL Immunology DB + CATNAP), 100% local.

Reimplemented (not imported) from the standalone project's
``lanl_catnap_engine.py``. Pure pandas/csv logic, no subprocess: this
plugin wraps no external tool, only two local reference databases.

Does NOT do structural alignment or HXB2 coordinate mapping: candidate
linear sequences are compared directly against ``ab_all.csv``'s reported
epitopes by substring overlap. Deliberately simple (no network, no new
dependencies, pure pandas over already-downloaded CSVs) at the cost of not
capturing conformational epitopes (out of scope: would need 3D structure,
not just sequence).
"""

import csv
import re
from pathlib import Path
from typing import List, Optional

import pandas as pd

_OUTPUT_COLUMNS = [
    'sequence', 'antibody_name', 'epitope_sequence', 'match_length', 'epitope_name',
    'hxb2_location', 'neutralizing', 'antibody_type', 'subtype', 'binding_region',
    'catnap_mean_ic50', 'catnap_n_viruses',
]

_AA_ONLY = re.compile(r'[A-Zx]+')

_REQUIRED_LANL_COLUMNS = [
    'Antibody name (alias)', 'Epitope', 'Epitope name', 'HXB2 protein location',
    'Neutralizing', 'Antibody type', 'Binding region', 'Subtype',
]

_REQUIRED_CATNAP_COLUMNS = ['Name', 'Mean panel IC50', '# of viruses tested']


class LANLCATNAPParseError(Exception):
    """The LANL/CATNAP reference file does not match the expected format."""


def load_bnab_epitopes(lanlAbAllPath: Path) -> pd.DataFrame:
    """Parse ``ab_all.csv``, keeping only records with a linear epitope (single-chain AA sequence)."""
    with open(lanlAbAllPath, encoding='utf-8', errors='replace', newline='') as fh:
        reader = csv.reader(fh)
        header = next(reader)
        rows = list(reader)

    idx = {name: i for i, name in enumerate(header)}
    missing = [c for c in _REQUIRED_LANL_COLUMNS if c not in idx]
    if missing:
        raise LANLCATNAPParseError(f"'{lanlAbAllPath}' does not have the expected columns: missing {missing}.")

    records = []
    for row in rows:
        epitope = row[idx['Epitope']].strip()
        if not epitope or not _AA_ONLY.fullmatch(epitope):
            continue  # conformational epitope, composite notation ('A + B'), or empty -- out of scope
        records.append({
            'antibody_name': row[idx['Antibody name (alias)']].strip(),
            'epitope_sequence': epitope.upper(),
            'epitope_name': row[idx['Epitope name']].strip(),
            'hxb2_location': row[idx['HXB2 protein location']].strip(),
            'neutralizing': row[idx['Neutralizing']].strip(),
            'antibody_type': row[idx['Antibody type']].strip(),
            'binding_region': row[idx['Binding region']].strip(),
            'subtype': row[idx['Subtype']].strip(),
        })
    return pd.DataFrame.from_records(records)


def load_catnap_potency(catnapAbsPath: Path) -> pd.DataFrame:
    """Parse ``abs_YYYY-MM-DD.txt`` (CATNAP) to append neutralization potency/breadth per antibody."""
    raw = pd.read_csv(catnapAbsPath, sep='\t', dtype=str)
    missing = [c for c in _REQUIRED_CATNAP_COLUMNS if c not in raw.columns]
    if missing:
        raise LANLCATNAPParseError(f"'{catnapAbsPath}' does not have the expected columns: missing {missing}.")

    result = raw[_REQUIRED_CATNAP_COLUMNS].copy()
    result.columns = ['antibody_name_norm', 'catnap_mean_ic50', 'catnap_n_viruses']
    result['antibody_name_norm'] = result['antibody_name_norm'].str.strip().str.upper()
    result['catnap_mean_ic50'] = pd.to_numeric(result['catnap_mean_ic50'], errors='coerce')
    result['catnap_n_viruses'] = pd.to_numeric(result['catnap_n_viruses'], errors='coerce')
    return result.drop_duplicates(subset='antibody_name_norm', keep='first')


def longest_common_substring_len(a: str, b: str) -> int:
    """Length of the longest common substring between ``a`` and ``b`` (DP O(len(a)*len(b)))."""
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        curr = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                curr[j] = prev[j - 1] + 1
                best = max(best, curr[j])
        prev = curr
    return best


def query_bnab_crossref(
    sequences: List[str], lanlAbAllPath: Path, catnapAbsPath: Optional[Path] = None,
    minOverlap: int = 6,
) -> pd.DataFrame:
    """Cross-reference ``sequences`` against known bnAb linear epitopes (LANL Immunology DB).

    Args:
        sequences: Candidate peptides/sequences to evaluate.
        lanlAbAllPath: Path to LANL's ``ab_all.csv``.
        catnapAbsPath: Path to CATNAP's ``abs_YYYY-MM-DD.txt`` (optional).
        minOverlap: Minimum substring overlap length to report a match. For
            reference epitopes SHORTER than this, the full epitope match is
            required instead (never a looser threshold than the epitope
            itself).

    Returns:
        DataFrame with one row per (candidate, reference epitope) pair that
        overlaps enough, columns ``_OUTPUT_COLUMNS``. Empty if no match or
        ``sequences`` is empty.
    """
    if not sequences:
        return pd.DataFrame(columns=_OUTPUT_COLUMNS)

    bnabDf = load_bnab_epitopes(lanlAbAllPath)
    if bnabDf.empty:
        return pd.DataFrame(columns=_OUTPUT_COLUMNS)

    potencyDf = load_catnap_potency(catnapAbsPath) if catnapAbsPath is not None else None

    rows = []
    for seq in sequences:
        seqUpper = seq.upper()
        for ref in bnabDf.itertuples(index=False):
            requiredOverlap = min(minOverlap, len(ref.epitope_sequence))
            matchLen = longest_common_substring_len(seqUpper, ref.epitope_sequence)
            if matchLen < requiredOverlap:
                continue

            record = {
                'sequence': seq,
                'antibody_name': ref.antibody_name,
                'epitope_sequence': ref.epitope_sequence,
                'match_length': matchLen,
                'epitope_name': ref.epitope_name,
                'hxb2_location': ref.hxb2_location,
                'neutralizing': ref.neutralizing,
                'antibody_type': ref.antibody_type,
                'subtype': ref.subtype,
                'binding_region': ref.binding_region,
                'catnap_mean_ic50': pd.NA,
                'catnap_n_viruses': pd.NA,
            }
            if potencyDf is not None:
                hit = potencyDf[potencyDf['antibody_name_norm'] == ref.antibody_name.strip().upper()]
                if not hit.empty:
                    record['catnap_mean_ic50'] = hit.iloc[0]['catnap_mean_ic50']
                    record['catnap_n_viruses'] = hit.iloc[0]['catnap_n_viruses']
            rows.append(record)

    return pd.DataFrame(rows, columns=_OUTPUT_COLUMNS) if rows else pd.DataFrame(columns=_OUTPUT_COLUMNS)
