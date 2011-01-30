# ----------------------------------------------------
# File: engrid-vtk-win_paraview.pri
# Date: 29 January 2011
# Author: Philippose Rajan
#
# This file contains the configuration required to 
# compile Engrid using the VTK version which is part 
# of ParaView compiled from sources on the Windows 
# platform using Microsoft Visual C++ (Express Edition).
#
# Note #1: This requires the option "INSTALL_DEVELOPMENT"
# to be enabled in the ParaView CMake configuration. 
#
# NOTE #2: The location of the ParaView binary installation 
# and the Paraview version need to be  specified as environment 
# variables:
# PV_INSTDIR
# PV_VER (Ex. 3.9)
# -----------------------------------------------------

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

