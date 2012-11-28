// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
// +                                                                      +
// + enGrid is free software: you can redistribute it and/or modify       +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
#ifndef ENGRID_VERSION
#if defined(WIN64) || defined(__MAC64) || defined(__LINUX64)
  #define ENGRID_VERSION "1.4.0_x64"
#else
  #define ENGRID_VERSION "1.4.0"
#endif
#define ENGRID_VERSION_FULLVER 1,4,0,0
#define ENGRID_COMPANY_NAME "enGits GmbH"
#define ENGRID_DESCRIPTION "enGrid is an open-source mesh generation software with CFD applications in mind."
#define ENGRID_COPYRIGHT "GNU Public License (GPL). Developed by enGits GmbH - http://engits.eu"
#define ENGRID_NAME "enGrid"
#endif
