# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2014 enGits GmbH                                      +
# +                                                                      +
# + enGrid is free software: you can redistribute it and/or modify       +
# + it under the terms of the GNU General Public License as published by +
# + the Free Software Foundation, either version 3 of the License, or    +
# + (at your option) any later version.                                  +
# +                                                                      +
# + enGrid is distributed in the hope that it will be useful,            +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
# + GNU General Public License for more details.                         +
# +                                                                      +
# + You should have received a copy of the GNU General Public License    +
# + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Default values for the above variables
DEF_PV_INSTDIR = "I:/ParaView_Git/ParaView-bin"
DEF_PV_VERSION = "3.11"

# ------ Check and assign the right values to the variables ------
TMPVAR = $$(PV_VER)
!isEmpty(TMPVAR) {
    PV_VERSION = $(PV_VER)
} else {
    PV_VERSION = $${DEF_PV_VERSION}
}

TMPVAR = $$(PV_INSTDIR)
!isEmpty(TMPVAR) {
    VTK_LIB_DIR    = $(PV_INSTDIR)/lib/paraview-$${PV_VERSION}
    VTK_PV_INC_DIR = $(PV_INSTDIR)/include/paraview-$${PV_VERSION}
} else {
    VTK_LIB_DIR    = $${DEF_PV_INSTDIR}/lib/paraview-$${PV_VERSION}
    VTK_PV_INC_DIR = $${DEF_PV_INSTDIR}/include/paraview-$${PV_VERSION}    
}
# -----------------------------------------------------------------

INCLUDEPATH += $${VTK_PV_INC_DIR}

LIBS += -L$${VTK_LIB_DIR}

