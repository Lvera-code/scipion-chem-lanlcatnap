# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Enzo Sierra (enzogael57@gmail.com)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This protocol is used to cross-reference peptide candidates against known
HIV broadly-neutralizing-antibody (bnAb) linear epitopes, using the local
LANL Immunology DB (+ optional CATNAP neutralization data).
"""

import os

import pandas as pd
from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Boolean, Integer, String
from pyworkflow.protocol import params

from .. import Plugin as lanlcatnapPlugin
from ..constants import DEFAULT_MIN_OVERLAP, LANLCATNAP_DIC
from ..utils.lanlcatnap import query_bnab_crossref


class ProtLANLCATNAPCrossref(EMProtocol):
    """
    AI Generated:

    Cross-references every input ROI's peptide against known HIV
    broadly-neutralizing-antibody (bnAb) linear epitopes (LANL Immunology
    DB, optionally enriched with CATNAP neutralization potency), and
    annotates (does NOT filter) each ROI with a summary of its best match.
    Purely informative for any pipeline: only relevant when the input is
    HIV Env, but harmless (zero matches) otherwise.

    Does NOT do structural alignment or HXB2 coordinate mapping: compares
    candidate linear sequences directly against LANL's reported epitopes
    by longest-common-substring overlap. Wraps no external tool: pure
    pandas/csv logic over two local reference databases.

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, annotated with
    ``_bnabMatchCount`` (int, total matching reference epitopes),
    ``_bnabNeutralizingMatch`` (bool, True if >=1 match is a confirmed
    neutralizing antibody), ``_bnabBestAntibody``/``_bnabBestMatchLength``
    (the single longest match, name and length). The full detail (every
    match, not just the best one) is persisted to
    ``extra/bnab_crossref.csv``.
    """

    _label = 'lanl-catnap bnab crossref'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Peptide candidates to cross-reference against known bnAb epitopes.')
        form.addParam('minOverlap', params.IntParam, default=DEFAULT_MIN_OVERLAP,
                       label='Min. substring overlap (aa): ',
                       help='Minimum substring overlap to report a match. Reference epitopes '
                            'shorter than this still require their own full length as the match.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.crossrefStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getCrossrefPath(self):
        return self._getExtraPath('bnab_crossref.csv')

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def crossrefStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        catnapPath = lanlcatnapPlugin.getVar(LANLCATNAP_DIC['catnap_abs_path']) or None
        crossrefDf = query_bnab_crossref(
            sequences, lanlAbAllPath=lanlcatnapPlugin.getVar(LANLCATNAP_DIC['ab_all_path']),
            catnapAbsPath=catnapPath, minOverlap=self.minOverlap.get(),
        )
        crossrefDf.to_csv(self._getCrossrefPath(), index=False)

    def createOutputStep(self):
        rois = self._getRois()
        crossrefDf = pd.read_csv(self._getCrossrefPath()) if os.path.isfile(self._getCrossrefPath()) else pd.DataFrame()

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi in rois:
            matches = crossrefDf[crossrefDf['sequence'] == roi.getROISequence()] if not crossrefDf.empty else crossrefDf

            roi._bnabMatchCount = Integer(len(matches))
            roi._bnabNeutralizingMatch = Boolean(bool((matches['neutralizing'] == 'yes').any()) if len(matches) else False)
            if len(matches):
                best = matches.sort_values('match_length', ascending=False).iloc[0]
                roi._bnabBestAntibody = String(best['antibody_name'])
                roi._bnabBestMatchLength = Integer(int(best['match_length']))
            else:
                roi._bnabBestAntibody = String('')
                roi._bnabBestMatchLength = Integer(0)
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return lanlcatnapPlugin.validateInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                nMatch = sum(1 for roi in outROIs if roi._bnabMatchCount.get() > 0)
                nNeutralizing = sum(1 for roi in outROIs if roi._bnabNeutralizingMatch.get())
                summary.append(f'{nMatch}/{len(outROIs)} candidate(s) match >= 1 known bnAb epitope '
                               f'({nNeutralizing} with a confirmed neutralizing antibody).')
        return summary
