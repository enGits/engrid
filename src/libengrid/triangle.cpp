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
#include "triangle.h"
#include "geometrytools.h"
#include "engrid.h"
#include "utilities.h"
#include "egvtkobject.h"

Triangle::Triangle() : EgVtkObject()
{
  setDefaults();
  setupTriangle();
}

Triangle::Triangle(vec3_t a, vec3_t b, vec3_t c) : EgVtkObject()
{
  setDefaults();
  m_Xa = a;
  m_Xb = b;
  m_Xc = c;
  setupTriangle();
}

Triangle::Triangle(vtkUnstructuredGrid* grid, vtkIdType id_a, vtkIdType id_b, vtkIdType id_c) : EgVtkObject()
{
  setDefaults();
  m_IdA = id_a;
  m_IdB = id_b;
  m_IdC = id_c;
  grid->GetPoints()->GetPoint(id_a, m_Xa.data());
  grid->GetPoints()->GetPoint(id_b, m_Xb.data());
  grid->GetPoints()->GetPoint(id_c, m_Xc.data());
  setupTriangle();
}

Triangle::Triangle(vtkUnstructuredGrid* grid, vtkIdType id_cell)  : EgVtkObject()
{
  setDefaults();
  vtkIdType Npts, *pts;
  grid->GetCellPoints(id_cell, Npts, pts);
  if (Npts == 3) {
    m_IdA = pts[0];
    m_IdB = pts[1];
    m_IdC = pts[2];
    grid->GetPoints()->GetPoint(m_IdA, m_Xa.data());
    grid->GetPoints()->GetPoint(m_IdB, m_Xb.data());
    grid->GetPoints()->GetPoint(m_IdC, m_Xc.data());
    setupTriangle();
  } else {
    EG_ERR_RETURN("only triangles allowed at the moment");
  }
}

void Triangle::setDefaults()
{
  m_IdA = 0;
  m_IdB = 0;
  m_IdC = 0;
  m_Xa = vec3_t(0, 0, 0);
  m_Xb = vec3_t(0, 1, 0);
  m_Xc = vec3_t(0, 0, 1);
  m_NormalA = vec3_t(0, 0, 0);
  m_NormalB = vec3_t(0, 0, 0);
  m_NormalC = vec3_t(0, 0, 0);
}

void Triangle::setupTriangle()
{
  m_HasNeighbour.fill(false, 6);

  m_G1 = m_Xb - m_Xa;
  m_G2 = m_Xc - m_Xa;
  m_G3 = m_G1.cross(m_G2);
  
  if(m_G3.abs2() <= 0) {
    m_Valid = false;
  } else {
    m_Valid = true;
  }
  
  if(!checkVector(m_G3)) {
    qWarning() << "m_G1=" << m_G1;
    qWarning() << "m_G2=" << m_G2;
    qWarning() << "m_G3=" << m_G3;
    EG_BUG;
  }
  
  m_A  = 0.5 * m_G3.abs();
  if (m_Valid) {
    m_G3.normalise();
  }

  if(!checkVector(m_G3)) {
    qWarning() << "m_G1="<<m_G1;
    qWarning() << "m_G2="<<m_G2;
    qWarning() << "m_G3="<<m_G3;
    qWarning() << "m_Xa="<<m_Xa;
    qWarning() << "m_Xb="<<m_Xb;
    qWarning() << "m_Xc="<<m_Xc;
    this->saveTriangle("crash");
    EG_BUG;
  }
  
  m_G.column(0, m_G1);
  m_G.column(1, m_G2);
  m_G.column(2, m_G3);
  m_GI = m_G.inverse();

  m_SmallestLength = (m_Xb - m_Xa).abs();
  m_SmallestLength = min(m_SmallestLength, (m_Xc - m_Xb).abs());
  m_SmallestLength = min(m_SmallestLength, (m_Xa - m_Xc).abs());

  // compute minimal height
  double ha = (GeometryTools::projectPointOnEdge(m_Xa, m_Xb, (m_Xc - m_Xb)) - m_Xa).abs();
  double hb = (GeometryTools::projectPointOnEdge(m_Xb, m_Xa, (m_Xc - m_Xa)) - m_Xb).abs();
  double hc = (GeometryTools::projectPointOnEdge(m_Xc, m_Xa, (m_Xb - m_Xa)) - m_Xc).abs();
  m_SmallestHeight = min(ha, min(hb, hc));

}

vec3_t Triangle::local3DToGlobal3D(vec3_t r)
{
  return m_Xa + m_G*r;
}

vec3_t Triangle::global3DToLocal3D(vec3_t x)
{
  vec3_t tmp = x - m_Xa;
  return m_GI*tmp;
}

vec3_t Triangle::local2DToGlobal3D(vec2_t r)
{
  return local3DToGlobal3D(vec3_t(r[0], r[1], 0));
}

vec2_t Triangle::global3DToLocal2D(vec3_t x)
{
  vec3_t r = global3DToLocal3D(x);
  return vec2_t(r[0], r[1]);
}

bool Triangle::snapOntoTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d, int& side, bool restrict_to_triangle)
{
  side = -1;
  double scal = (xp - this->m_Xa) * this->m_G3;
  vec3_t x1, x2;
  if (scal > 0) {
    x1 = xp + this->m_G3;
    x2 = xp - scal * this->m_G3 - this->m_G3;
  } else {
    x1 = xp - this->m_G3;
    x2 = xp - scal * this->m_G3 + this->m_G3;
  }
  // (xi,ri) gets set to the intersection of the line with the plane here!
  bool intersects_face = GeometryTools::intersectEdgeAndTriangle(this->m_Xa, this->m_Xb, this->m_Xc, x1, x2, xi, ri);
  vec3_t xi_free = xi;
  if (intersects_face) {
    vec3_t dx = xp - this->m_Xa;
    d = fabs(dx * this->m_G3);
  } else {
    double kab = GeometryTools::intersection(this->m_Xa, this->m_Xb - this->m_Xa, xp, this->m_Xb - this->m_Xa);
    double kac = GeometryTools::intersection(this->m_Xa, this->m_Xc - this->m_Xa, xp, this->m_Xc - this->m_Xa);
    double kbc = GeometryTools::intersection(this->m_Xb, this->m_Xc - this->m_Xb, xp, this->m_Xc - this->m_Xb);

    double dab = (this->m_Xa + kab * (this->m_Xb - this->m_Xa) - xp).abs2();
    double dac = (this->m_Xa + kac * (this->m_Xc - this->m_Xa) - xp).abs2();
    double dbc = (this->m_Xb + kbc * (this->m_Xc - this->m_Xb) - xp).abs2();
    double da = (this->m_Xa - xp).abs2();
    double db = (this->m_Xb - xp).abs2();
    double dc = (this->m_Xc - xp).abs2();

    bool set = false;
    d = 1e99;//max(max(max(max(max(dab,dac),dbc),da),db),dc);

    if (dab < d) {
      if ((kab >= 0) && (kab <= 1)) {
        xi = this->m_Xa + kab * (this->m_Xb - this->m_Xa);
        ri = vec3_t(kab, 0, 0);
        d = dab;
        set = true;
        side = 0;
      }
    }
    if (dbc < d) {
      if ((kbc >= 0) && (kbc <= 1)) {
        xi = this->m_Xb + kbc * (this->m_Xc - this->m_Xb);
        ri = vec3_t(1 - kbc, kbc, 0);
        d = dbc;
        set = true;
        side = 1;
      }
    }
    if (dac < d) {
      if ((kac >= 0) && (kac <= 1)) {
        xi = this->m_Xa + kac * (this->m_Xc - this->m_Xa);
        ri = vec3_t(0, kac, 0);
        d = dac;
        set = true;
        side = 2;
      }
    }
    if (da < d) {
      xi = this->m_Xa;
      ri = vec3_t(0, 0);
      d = da;
      set = true;
      side = 3;
    }
    if (db < d) {
      xi = this->m_Xb;
      ri = vec3_t(1, 0);
      d = db;
      set = true;
      side = 4;
    }
    if (dc < d) {
      xi = this->m_Xc;
      ri = vec3_t(0, 1);
      d = dc;
      set = true;
      side = 5;
    }
    if (!set) {
      EG_BUG;
    }
    d = sqrt(d);
  }
  if (!intersects_face && !restrict_to_triangle) {
    xi = xi_free;
  }
  if (xi[0] > 1e98) { // should never happen
    EG_BUG;
  }
  return intersects_face;
}

void Triangle::saveTriangle(QString filename)
{
  int N_cells = 1;
  int N_points = 3;
  
  EG_VTKSP(vtkUnstructuredGrid, triangle_grid);
  allocateGrid(triangle_grid , N_cells, N_points);
  
  vtkIdType node_count = 0;
  int cell_count = 0;
  
  vtkIdType pts[3];
  triangle_grid->GetPoints()->SetPoint(node_count, m_Xa.data()); pts[0]=node_count; node_count++;
  triangle_grid->GetPoints()->SetPoint(node_count, m_Xb.data()); pts[1]=node_count; node_count++;
  triangle_grid->GetPoints()->SetPoint(node_count, m_Xc.data()); pts[2]=node_count; node_count++;
  
  triangle_grid->InsertNextCell(VTK_TRIANGLE,3,pts);cell_count++;
  
  saveGrid(triangle_grid, filename+"_triangle_grid");
}

void Triangle::setNormals(vec3_t na, vec3_t nb, vec3_t nc)
{
  m_NormalA = na;
  m_NormalB = nb;
  m_NormalC = nc;
  // compute normal vectors in local coordinate system
  m_RNormalA = global3DToLocal3D(a() + nA());
  m_RNormalB = global3DToLocal3D(a() + nB());
  m_RNormalC = global3DToLocal3D(a() + nC());
}
