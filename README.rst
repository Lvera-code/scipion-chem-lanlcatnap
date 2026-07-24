================================
LANL/CATNAP bnAb cross-reference Scipion plugin
================================

Scipion framework plugin cross-referencing peptide candidates against
known HIV broadly-neutralizing-antibody (bnAb) linear epitopes, using the
local LANL HIV Molecular Immunology Database (+ optional CATNAP
neutralization potency data, Yoon et al. 2015).

The plugin implements a single protocol, ``ProtLANLCATNAPCrossref``, which
annotates (does **not** filter) every input ROI with a summary of its best
matching bnAb epitope (``_bnabMatchCount``/``_bnabNeutralizingMatch``/
``_bnabBestAntibody``/``_bnabBestMatchLength``). Purely informative for
any pipeline: only relevant when the input is HIV Env, but harmless (zero
matches) otherwise. The full match detail (every match, not just the
single best one) is persisted to ``extra/bnab_crossref.csv``.

This plugin wraps **no external tool**: pure pandas/csv logic over two
local reference databases, not bundled with the plugin (LANL/CATNAP's own
terms of use do not clearly permit redistribution -- download once,
reference by path in ``scipion.conf``, same convention as every other
academic dataset/tool in this project).

================================
Manual setup
================================

1. Download the LANL HIV Molecular Immunology DB's ``ab_all.csv`` from
   https://www.hiv.lanl.gov/content/immunology/ (required).
2. Optionally, download CATNAP's ``abs_YYYY-MM-DD.txt`` from
   https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/ (used
   only to append neutralization potency/breadth when available).
3. In ``scipion.conf``, set:

.. code-block::

      LANL_AB_ALL_PATH = <path to ab_all.csv>
      CATNAP_ABS_PATH = <path to abs_YYYY-MM-DD.txt>   # optional

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-lanlcatnap.git
            cd scipion-chem-lanlcatnap
            scipion3 installp -p . --devel
