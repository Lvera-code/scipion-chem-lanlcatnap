=========
CHANGES
=========

0.2.0
=====
- Deprecated: this plugin has been folded into ``pwchem`` core, following
  the same treatment previously applied to EpiDope. ``ProtLANLCATNAPCrossref``
  now lives in ``pwchem/protocols/Sequences/protocol_lanlcatnap_crossref.py``.
  See ``README.rst`` for details. This repository is kept for reference
  only and should not be installed.

0.1.0
=====
- Initial release: LANL Immunology DB + CATNAP bnAb cross-reference
  annotation protocol (``ProtLANLCATNAPCrossref``). Wraps no external
  tool: pure pandas/csv logic over two local reference databases.
