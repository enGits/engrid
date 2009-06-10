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
  
  bool    m_Clip;
  bool    m_ExtrTetras;
  bool    m_ExtrHexes;
  bool    m_ExtrWedges;
  bool    m_ExtrPyramids;
  vec3_t  m_X;
  vec3_t  m_N;

public: // methods
  
  static vtkEgExtractVolumeCells* New();
  void SetX(vec3_t x);
  void Setx(double x);
  void Sety(double y);
  void Setz(double z);
  void SetN(vec3_t n);
  void Setnx(double nx);
  void Setny(double ny);
  void Setnz(double nz);
  void SetClippingOn();
  void SetClippingOff();
  void SetAllOn();
  void SetAllOff();
  void SetTetrasOn();
  void SetTetrasOff();
  void SetPyramidsOn();
  void SetPyramidsOff();
  void SetWedgesOn();
  void SetWedgesOff();
  void SetHexesOn();
  void SetHexesOff();
  
protected: // methods
  
  vtkEgExtractVolumeCells();
  ~vtkEgExtractVolumeCells() {}
  virtual void ExecuteEg();
  
private: // methods
  
  vtkEgExtractVolumeCells (const vtkEgExtractVolumeCells&);
  void operator= (const vtkEgExtractVolumeCells&);
  
};

#endif
