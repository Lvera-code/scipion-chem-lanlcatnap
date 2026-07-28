================================
DEPRECATED: LANL/CATNAP bnAb cross-reference Scipion plugin
================================

This plugin is **deprecated** and pending archival. Following the same
reasoning applied earlier to ``scipion-chem-bcellepitope``, and review
feedback from Blanca (CNB Biocomputing unit), this protocol has been
folded into ``pwchem`` core instead of living in its own plugin: it wraps
no external tool/binary (pure pandas/csv logic over two local reference
databases), so it needs no dedicated conda environment and does not
justify a standalone plugin repo, same class of decision as EpiDope.

Where it moved to
================================

``ProtLANLCATNAPCrossref`` now lives in ``pwchem``
(``pwchem/protocols/Sequences/protocol_lanlcatnap_crossref.py``, same
place as EpiDope). Configuration is unchanged: set ``LANL_AB_ALL_PATH``
(required) and ``CATNAP_ABS_PATH`` (optional) in ``scipion.conf`` -- see
that plugin's own docs for the manual download steps (LANL HIV Molecular
Immunology DB + optional CATNAP neutralization data). This stays a manual
download, unlike most other tools recently migrated to auto-install:
verified 2026-07-28 that hiv.lanl.gov has no stable download URL for its
antibody DB (only an interactive search form) and explicitly states it is
increasing its blocking of automated traffic.

Status
================================

Not yet PR'd upstream to ``scipion-chem``/pwchem's real GitHub org (same
status as the EpiDope migration this follows). Until then, this
repository is kept around for reference only and should not be installed.
