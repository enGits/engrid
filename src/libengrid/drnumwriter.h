// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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

#ifndef DRNUMWRITER_H
#define DRNUMWRITER_H

#include "iooperation.h"
#include "edgelengthsourcemanager.h"

#include <QFile>
#include <QTextStream>

/**
 * @brief A very experimental export function for DrNUM grids.
 */
class DrNumWriter : public IOOperation
{

protected: // data types

  struct cart_patch_t
  {
    vec3_t    x0, gi, gj;
    int       Ni, Nj, Nk, sl_i1, sl_i2, sl_j1, sl_j2, sl_k1, sl_k2;
    double    Li, Lj, Lk;
    QString   fx, fy, fz, bx, bX, by, bY, bz, bZ, s;
    vtkIdType id_cell;
  };


protected: // attributes

  QString m_PatchFile;   ///< the name of the DrNUM patches file
  QString m_ComplexPath; ///< the path to the DrNUM directory for complex patches

  QVector<cart_patch_t>   m_CartPatches;
  QVector<int>            m_CellToCartPatch;
  int                     m_OverlapLayers;
  QVector<double>         m_H;
  double                  m_MaxEdgeLength;
  double                  m_MinEdgeLength;
  double                  m_GrowthFactor;
  EdgeLengthSourceManager m_ELSManager;



protected: // methods

  void computeMeshDensity();
  void createAllCartPatches();
  void writeCartPatches(QTextStream &s);
  QString boundaryCode(vtkIdType id_cell, int i);

  virtual void operate();


public: // methods

  DrNumWriter();

};

#endif // DRNUMWRITER_H
