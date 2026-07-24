from pwem.protocols import ProtImportSequence
from pwchem.protocols import ProtDefineSeqROI
from pyworkflow.tests import BaseTest, setupTestProject

from ..protocols import ProtLANLCATNAPCrossref


class TestLANLCATNAPCrossref(BaseTest):
    NAME = 'LANLCATNAP_TEST_SEQ'
    DESCRIPTION = 'GP120 MPER region + an unrelated N-term fragment'
    # PEP1: real HIV-1 Env gp41 MPER region (matches this project's GP120
    # fixture sequence used across other test suites). PEP2: an unrelated
    # real GP120 N-terminal fragment, included as a genuine negative
    # control (0 real matches against LANL's bnAb epitopes).
    PEPTIDES = ['WASLWNWFNITNWLWYIKIFIMIVGGLVGLRIVFAVLSIVNRV', 'RVKEKYQHL']
    SPACER = 'GGG'
    AMINOACIDSSEQ = SPACER.join(PEPTIDES)

    # Real output of query_bnab_crossref against the real local
    # reference_db/lanl_immunology/ab_all.csv + reference_db/catnap/
    # abs_2026-07-01.txt (not estimated): PEP1 matches 9 known bnAb
    # epitopes (DH416/DH415/DH413/Z13e1/WR316/VRC46.01/VRC43.01/10E8v4/
    # 10E8), several confirmed neutralizing (Z13e1, VRC46.01, VRC43.01,
    # 10E8v4, 10E8); the single LONGEST match is Z13e1 (8aa, real
    # CATNAP-confirmed neutralizing antibody). PEP2 matches nothing.
    EXPECTED = {
        'WASLWNWFNITNWLWYIKIFIMIVGGLVGLRIVFAVLSIVNRV': {
            'matchCount': 9, 'neutralizing': True, 'bestAntibody': 'Z13e1', 'bestMatchLength': 8,
        },
        'RVKEKYQHL': {
            'matchCount': 0, 'neutralizing': False, 'bestAntibody': '', 'bestMatchLength': 0,
        },
    }

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls._runImportSeq()
        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)

        cls.protSeedROIs = cls._runDefSeqROIs(cls.protImportSeq)
        cls._waitOutput(cls.protSeedROIs, 'outputROIs', sleepTime=5)

    @classmethod
    def _runImportSeq(cls):
        kwargs = {
            'inputSequenceName': cls.NAME,
            'inputSequenceDescription': cls.DESCRIPTION,
            'inputRawSequence': cls.AMINOACIDSSEQ,
        }
        cls.protImportSeq = cls.newProtocol(ProtImportSequence, **kwargs)
        cls.proj.launchProtocol(cls.protImportSeq, wait=False)

    @classmethod
    def _getWindows(cls):
        windows = []
        cursor = 0
        for pep in cls.PEPTIDES:
            start = cls.AMINOACIDSSEQ.index(pep, cursor) + 1
            end = start + len(pep) - 1
            windows.append((start, end))
            cursor = end
        return windows

    @classmethod
    def _runDefSeqROIs(cls, inProt):
        windows = cls._getWindows()
        inROIs = '\n'.join(
            '{}) Residues: {{"index": "{}-{}", "residues": "{}", "desc": "None"}}'.format(
                i, start, end, cls.AMINOACIDSSEQ[start - 1:end]
            )
            for i, (start, end) in enumerate(windows, 1)
        )
        protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
        protDefSeqROIs.inputSequence.set(inProt)
        protDefSeqROIs.inputSequence.setExtended('outputSequence')

        cls.proj.launchProtocol(protDefSeqROIs, wait=False)
        return protDefSeqROIs

    def test(self):
        protCrossref = self.newProtocol(ProtLANLCATNAPCrossref)
        protCrossref.inputROIs.set(self.protSeedROIs)
        protCrossref.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protCrossref, wait=True)

        outROIs = getattr(protCrossref, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertEqual(len(outROIs), len(self.PEPTIDES))

        for roi in outROIs:
            seq = roi.getROISequence()
            expected = self.EXPECTED[seq]
            self.assertEqual(roi._bnabMatchCount.get(), expected['matchCount'])
            self.assertEqual(roi._bnabNeutralizingMatch.get(), expected['neutralizing'])
            self.assertEqual(roi._bnabBestAntibody.get(), expected['bestAntibody'])
            self.assertEqual(roi._bnabBestMatchLength.get(), expected['bestMatchLength'])
