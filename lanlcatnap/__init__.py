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
This package contains a protocol for cross-referencing peptide candidates
against known HIV broadly-neutralizing-antibody (bnAb) linear epitopes,
using the local LANL Immunology DB (+ optional CATNAP neutralization data).
"""

import os

from pwchem import Plugin as pwchemPlugin

from .constants import LANLCATNAP_DIC, NOINSTALL_WARNING

_references = ['Yoon2015']


class Plugin(pwchemPlugin):
    """This plugin wraps no external tool/binary: pure pandas logic over
    two local reference databases, never auto-downloaded (LANL/CATNAP's
    own terms of use do not clearly permit redistribution). See
    ``validateInstallation`` and ``README.rst`` for the manual setup."""

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(LANLCATNAP_DIC['ab_all_path'], '')
        cls._defineVar(LANLCATNAP_DIC['catnap_abs_path'], '')

    @classmethod
    def defineBinaries(cls, env):
        """No-op: see class docstring."""
        pass

    @classmethod
    def validateInstallation(cls):
        """Check that this plugin's requirements are met. Returns a list of
        actionable error messages, empty if the installation is correct.
        CATNAP is optional (checked only if configured): only the LANL
        ab_all.csv path is required."""
        errors = []

        abAllPath = cls.getVar(LANLCATNAP_DIC['ab_all_path'])
        if not abAllPath or not os.path.isfile(abAllPath):
            errors.append(f"LANL_AB_ALL_PATH is not set or does not exist: '{abAllPath}'.")

        catnapPath = cls.getVar(LANLCATNAP_DIC['catnap_abs_path'])
        if catnapPath and not os.path.isfile(catnapPath):
            errors.append(f"CATNAP_ABS_PATH is set but does not exist: '{catnapPath}'.")

        if errors:
            errors.append(NOINSTALL_WARNING)
        return errors
