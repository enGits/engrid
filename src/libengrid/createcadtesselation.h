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
#ifndef CREATECADTESSELATION_H
#define CREATECADTESSELATION_H

#include "operation.h"
#include "cadinterface.h"

class CreateCadTesselation : public Operation
{

private: // attributes

  vec3_t m_X1;
  vec3_t m_X2;
  vec3_t m_XScan1;
  vec3_t m_XScan2;
  int    m_Ni;
  int    m_Nj;
  int    m_Nk;
  double m_Dx;
  double m_Dy;
  double m_Dz;
  double m_ScanMemory;
  bool   m_GeometryFound;
  int    m_NumIterations;
  int    m_PreservationType;
  double m_SmallestFeatureSize;
  double m_SmallestResolution;
  double m_TargetReduction;

  CadInterface* m_CadInterface;

protected: // methods

  virtual void operate();

  bool shootRay(vec3_t x, vec3_t v, vec3_t &x_in, vec3_t &x_out, vec3_t &n_in, vec3_t &n_out);
  void scan(bool create_grid, int interlaces = 0);

  double getx(int i) { return m_X1[0] + i*m_Dx; }
  double gety(int j) { return m_X1[1] + j*m_Dy; }
  double getz(int k) { return m_X1[2] + k*m_Dz; }
  vec3_t getX(int i, int j, int k) { return vec3_t(getx(i), gety(j), getz(k)); }
  int    getIdx(int i, int j, int k) { return i + j*m_Ni + k*m_Ni*m_Nj; }


public: // methods

  CreateCadTesselation(CadInterface* cad_interface);

  void setScanMemory(double mem) { m_ScanMemory = mem; }
  void setPreservationOff() { m_PreservationType = 0; }
  void setSolidPreservation() { m_PreservationType = 1; }
  void setFluidPreservation() { m_PreservationType = 2; }
  void setPreservationType(int t) { m_PreservationType = t; }
  void setSmoothingIterations(int n) { m_NumIterations = n; }
  bool preserveSolid() { return m_PreservationType == 1; }
  bool preserveFluid() { return m_PreservationType == 2; }
  void setSmallestFeatureSize(double sfs) { m_SmallestFeatureSize = sfs; }
  void setSmallestResolution(double h) { m_SmallestResolution = h; }
  void setTargetReduction(double tr) { m_TargetReduction = tr; }

};

#endif // CREATECADTESSELATION_H
