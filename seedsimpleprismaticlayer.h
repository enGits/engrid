//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#ifndef seedsimpleprismaticlayer_H
#define seedsimpleprismaticlayer_H

class SeedSimplePrismaticLayer;

#include "operation.h"

/**
 * The computation of the vectors providing the initial move away from the wall is not 
 * necessarily always correct.
 * It shoud be weighted with the angle of adjacent edges and then normalised with an
 * average edge length.
 */
class SeedSimplePrismaticLayer : public Operation
{
  
protected: // attributes
  
  QVector<vtkIdType>           layer_cells;
  QList<vtkIdType>             remesh_tetras;
  QVector<vtkIdType>           vol_cells;
  QVector<vtkIdType>           prismatic_cells;
  QVector<QVector<vtkIdType> > faces;
  QVector<vtkIdType>           old2new;
  
  int    N_new_points;
  int    N_new_cells;
  int    new_layer;
  int    old_layer;
  double layer_g;
  double layer_dg;
  
public: // methods
  
  SeedSimplePrismaticLayer();
  void setLayerCells(const QVector<vtkIdType> &cells);
  void getLayerCells(QVector<vtkIdType> &cells);
  void increaseLayerG(double d) { layer_g += d; layer_dg = d; };
  
protected: // methods
  
  void prepareLayer();
  void createBoundaryElements(vtkUnstructuredGrid *new_grid);
  int countBoundaryElements();
    
  virtual void operate();
  
};

#endif
