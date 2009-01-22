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
#ifndef __vtkEgExtractVolumeCells_h
#define __vtkEgExtractVolumeCells_h

class vtkEgExtractVolumeCells;

#include "vtkEgGridFilter.h"

class vtkEgExtractVolumeCells : public vtkEgGridFilter
{
  
protected: // attributes
  
  bool Clip;
  vec3_t X;
  vec3_t N;
  bool tetra;
  bool hexa;
  bool wedge;
  bool pyramid;
  
public: // methods
  
  static vtkEgExtractVolumeCells* New();
  void SetX(vec3_t x)   { X = x; Modified(); };
  void SetN(vec3_t n)   { N = n; Modified(); };
  void SetClippingOn()  { Clip = true;  Modified(); };
  void SetClippingOff() { Clip = false; Modified(); };
  void SetAllOn()       { tetra = true; hexa = true; wedge = true; pyramid = true; };
  void SetAllOff()      { tetra = false; hexa = false; wedge = false; pyramid = false; };
  void SetTetrasOn()    { tetra = true; };
  void SetTetrasOff()   { tetra = false; };
  void SetPyramidsOn()  { pyramid = true; };
  void SetPyramidsOff() { pyramid = false; };
  void SetWedgesOn()    { wedge = true; };
  void SetWedgesOff()   { wedge = false; };
  void SetHexasOn()     { hexa = true; };
  void SetHexasOff()    { hexa = false; };
  
protected: // methods
  
  vtkEgExtractVolumeCells();
  ~vtkEgExtractVolumeCells() {};
  virtual void ExecuteEg();
  
private: // methods
  
  vtkEgExtractVolumeCells (const vtkEgExtractVolumeCells&);
  void operator= (const vtkEgExtractVolumeCells&);
  
};

#endif
