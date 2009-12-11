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
#include "triangle.h"
#include "geometrytools.h"
#include "engrid.h"
#include "utilities.h"
#include "egvtkobject.h"

Triangle::Triangle() : EgVtkObject() {
  setDefaults();
  setupTriangle();
}

Triangle::Triangle(vec3_t a_a, vec3_t a_b, vec3_t a_c) : EgVtkObject() {
  setDefaults();
  m_a = a_a;
  m_b = a_b;
  m_c = a_c;
  setupTriangle();
}

Triangle::Triangle(vtkUnstructuredGrid* a_grid, vtkIdType a_id_a, vtkIdType a_id_b, vtkIdType a_id_c) : EgVtkObject() {
  setDefaults();
  m_id_a = a_id_a;
  m_id_b = a_id_b;
  m_id_c = a_id_c;
  a_grid->GetPoints()->GetPoint(m_id_a, m_a.data());
  a_grid->GetPoints()->GetPoint(m_id_b, m_b.data());
  a_grid->GetPoints()->GetPoint(m_id_c, m_c.data());
  setupTriangle();
}

Triangle::Triangle(vtkUnstructuredGrid* a_grid, vtkIdType a_id_cell)  : EgVtkObject() {
  setDefaults();
  vtkIdType Npts, *pts;
  a_grid->GetCellPoints(a_id_cell, Npts, pts);
  if (Npts == 3) {
    m_id_a = pts[0];
    m_id_b = pts[1];
    m_id_c = pts[2];
    a_grid->GetPoints()->GetPoint(m_id_a, m_a.data());
    a_grid->GetPoints()->GetPoint(m_id_b, m_b.data());
    a_grid->GetPoints()->GetPoint(m_id_c, m_c.data());
    setupTriangle();
  } else {
    EG_ERR_RETURN("only triangles allowed at the moment");
  }
}

void Triangle::setDefaults()
{
  m_id_a = 0;
  m_id_b = 0;
  m_id_c = 0;
  m_a = vec3_t(0, 0, 0);
  m_b = vec3_t(0, 1, 0);
  m_c = vec3_t(0, 0, 1);
  m_Normal_a = vec3_t(0, 0, 0);
  m_Normal_b = vec3_t(0, 0, 0);
  m_Normal_c = vec3_t(0, 0, 0);
}

void Triangle::setupTriangle() {
  m_has_neighbour.resize(6);
  m_has_neighbour[0] = false;
  m_has_neighbour[1] = false;
  m_has_neighbour[2] = false;
  m_has_neighbour[3] = false;
  m_has_neighbour[4] = false;
  m_has_neighbour[5] = false;

  m_g1 = m_b - m_a;
  m_g2 = m_c - m_a;
  m_g3 = m_g1.cross(m_g2);
  
  if(m_g3.abs2()<=0) {
    m_Valid = false;
  }
  else {
    m_Valid = true;
  }
  
  if(!checkVector(m_g3)) {
    qWarning()<<"m_g1="<<m_g1;
    qWarning()<<"m_g2="<<m_g2;
    qWarning()<<"m_g3="<<m_g3;
    vec3_t foo = vec3_t(0,0,0);
    vec3_t bar = vec3_t(0,0,0);
    qWarning()<<"foo.cross(bar)="<<foo.cross(bar);
    EG_BUG;
  }
  
  m_A  = 0.5 * m_g3.abs();
  if(m_Valid) m_g3.normalise();

  if(!checkVector(m_g3)) {
    qWarning()<<"m_g1="<<m_g1;
    qWarning()<<"m_g2="<<m_g2;
    qWarning()<<"m_g3="<<m_g3;
    qWarning()<<"m_a="<<m_a;
    qWarning()<<"m_b="<<m_b;
    qWarning()<<"m_c="<<m_c;
    this->saveTriangle("crash");
    EG_BUG;
  }
  
  m_G.column(0, m_g1);
  m_G.column(1, m_g2);
  m_G.column(2, m_g3);
  m_GI = m_G.inverse();

  m_smallest_length = (m_b - m_a).abs();
  m_smallest_length = min(m_smallest_length, (m_c - m_b).abs());
  m_smallest_length = min(m_smallest_length, (m_a - m_c).abs());
}

vec3_t Triangle::local3DToGlobal3D(vec3_t l_M) {
  return m_a + m_G*l_M;
}

vec3_t Triangle::global3DToLocal3D(vec3_t g_M) {
  vec3_t tmp = g_M - m_a;
  return m_GI*tmp;
}

vec3_t Triangle::local2DToGlobal3D(vec2_t l_M) {
  return local3DToGlobal3D(vec3_t(l_M[0], l_M[1], 0));
}

vec2_t Triangle::global3DToLocal2D(vec3_t g_M) {
  vec3_t l_M = global3DToLocal3D(g_M);
  return vec2_t(l_M[0], l_M[1]);
}

bool Triangle::projectOnTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d, int& side, bool restrict_to_triangle) {
  side = -1;
  double scal = (xp - this->m_a) * this->m_g3;
  vec3_t x1, x2;
  if (scal > 0) {
    x1 = xp + this->m_g3;
    x2 = xp - scal * this->m_g3 - this->m_g3;
  } else {
    x1 = xp - this->m_g3;
    x2 = xp - scal * this->m_g3 + this->m_g3;
  }
  // (xi,ri) gets set to the intersection of the line with the plane here!
  bool intersects_face = GeometryTools::intersectEdgeAndTriangle(this->m_a, this->m_b, this->m_c, x1, x2, xi, ri);
  if (intersects_face || !restrict_to_triangle) {
    vec3_t dx = xp - this->m_a;
    d = fabs(dx * this->m_g3);
  } else {
    double kab = GeometryTools::intersection(this->m_a, this->m_b - this->m_a, xp, this->m_b - this->m_a);
    double kac = GeometryTools::intersection(this->m_a, this->m_c - this->m_a, xp, this->m_c - this->m_a);
    double kbc = GeometryTools::intersection(this->m_b, this->m_c - this->m_b, xp, this->m_c - this->m_b);

    double dab = (this->m_a + kab * (this->m_b - this->m_a) - xp).abs();
    double dac = (this->m_a + kac * (this->m_c - this->m_a) - xp).abs();
    double dbc = (this->m_b + kbc * (this->m_c - this->m_b) - xp).abs();
    double da = (this->m_a - xp).abs();
    double db = (this->m_b - xp).abs();
    double dc = (this->m_c - xp).abs();

    bool set = false;
    d = 1e99;//max(max(max(max(max(dab,dac),dbc),da),db),dc);

    if (dab < d) {
      if ((kab >= 0) && (kab <= 1)) {
        xi = this->m_a + kab * (this->m_b - this->m_a);
        ri = vec3_t(kab, 0, 0);
        d = dab;
        set = true;
        side = 0;
      }
    }
    if (dbc < d) {
      if ((kbc >= 0) && (kbc <= 1)) {
        xi = this->m_b + kbc * (this->m_c - this->m_b);
        ri = vec3_t(1 - kbc, kbc, 0);
        d = dbc;
        set = true;
        side = 1;
      }
    }
    if (dac < d) {
      if ((kac >= 0) && (kac <= 1)) {
        xi = this->m_a + kac * (this->m_c - this->m_a);
        ri = vec3_t(0, kac, 0);
        d = dac;
        set = true;
        side = 2;
      }
    }
    if (da < d) {
      xi = this->m_a;
      ri = vec3_t(0, 0);
      d = da;
      set = true;
      side = 3;
    }
    if (db < d) {
      xi = this->m_b;
      ri = vec3_t(1, 0);
      d = db;
      set = true;
      side = 4;
    }
    if (dc < d) {
      xi = this->m_c;
      ri = vec3_t(0, 1);
      d = dc;
      set = true;
      side = 5;
    }
    if (!set) {
      EG_BUG;
    }
  }
  if (xi[0] > 1e98) { // should never happen
    EG_BUG;
  }
  /*  if (not( 0<=ri[0] && ri[0]<=1 && 0<=ri[1] && ri[1]<=1 && ri[2]==0 )) {
      qWarning()<<"ri="<<ri;
      EG_BUG;
    }*/
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
  triangle_grid->GetPoints()->SetPoint(node_count, m_a.data()); pts[0]=node_count; node_count++;
  triangle_grid->GetPoints()->SetPoint(node_count, m_b.data()); pts[1]=node_count; node_count++;
  triangle_grid->GetPoints()->SetPoint(node_count, m_c.data()); pts[2]=node_count; node_count++;
  
  triangle_grid->InsertNextCell(VTK_TRIANGLE,3,pts);cell_count++;
  
  saveGrid(triangle_grid, filename+"_triangle_grid");
}
