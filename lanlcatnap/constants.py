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

DEFAULT_VERSION = '1.0'

# This plugin wraps NO external tool/binary: it is pure pandas logic over
# two local, already-downloaded reference databases (LANL HIV Molecular
# Immunology DB + CATNAP). No subprocess, no venv, no conda env -- just two
# config variables pointing at the CSV/TSV files, which are NOT bundled
# with this plugin (large files, and LANL/CATNAP's own terms of use do not
# clearly permit redistribution -- the user downloads them once, same
# convention as every academic-license tool in this project, even though
# this data itself has no license restriction beyond that).
LANLCATNAP_DIC = {
    'name': 'LANLCATNAP',
    'version': DEFAULT_VERSION,
    'ab_all_path': 'LANL_AB_ALL_PATH',
    'catnap_abs_path': 'CATNAP_ABS_PATH',
}

READ_URL = 'https://github.com/Lvera-code/scipion-chem-lanlcatnap'
LANL_URL = 'https://www.hiv.lanl.gov/content/immunology/'
CATNAP_URL = 'https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/'

NOINSTALL_WARNING = (
    'Installation could not be completed because the local LANL Immunology '
    f'DB (ab_all.csv, download from {LANL_URL}) has not been found. CATNAP '
    f'({CATNAP_URL}) is optional (only used to append neutralization potency '
    'when available). Please check the scipion-chem-lanlcatnap README file '
    f'for more details: {READ_URL}'
)

# Below this, a substring overlap is statistical noise (appears by chance
# in any protein): see the real distribution of epitope lengths in
# ab_all.csv (min 3, median 12, max 47). For reference epitopes already
# SHORTER than this, the full epitope match is still required (never a
# looser threshold than the epitope itself).
DEFAULT_MIN_OVERLAP = 6
