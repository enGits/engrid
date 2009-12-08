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
#include "surfaceprojection.h"

#include "beziertriangle.h"
#include <math.h>

#include <vtkUnstructuredGridWriter.h>

///\todo check orientation of triangle on which we project to avoid projecting on the wrong side in case of close opposite sides of a surface
///\todo Delete those grids somewhere
SurfaceProjection::SurfaceProjection() : SurfaceAlgorithm()
{
  m_BGrid = vtkUnstructuredGrid::New();
  
  m_InterpolationGrid = vtkUnstructuredGrid::New();
  m_BezierGrid = vtkUnstructuredGrid::New();
  
  m_Relax = 0.9;
  m_DistWeight = 1.0;
  m_DistExp = 1.0;
  m_DirWeight = 1.0;
  m_DirExp = 1.0;
  m_WeightOffset = 0.001;
  m_MinOTLength = 0.0;
  m_MaxIter = 10;
  m_ConvLimit = 0.1;
  m_RadiusFactor = 0.2;
  m_UseLevelSet = false;
  double max_cells;
  m_MaxOTCells = 2000;
  m_NumDirect = 0;
  m_NumFull = 0;
  
  this->setGrid(m_BGrid);
  
  getSet("surface meshing", "correct curvature (experimental)", false, m_correctCurvature);
  
  m_ExactMode = 0;
}

SurfaceProjection::~SurfaceProjection()
{
//   m_InterpolationGrid->Delete();
//   m_BezierGrid->Delete();
}

void SurfaceProjection::setBackgroundGrid_initOctree()
{
  writeGrid(m_BGrid, "background");
  double bounds[6];
  m_BGrid->GetBounds(bounds);
  vec3_t x1(bounds[0], bounds[2], bounds[4]);
  vec3_t x2(bounds[1], bounds[3], bounds[5]);
  vec3_t xm = 0.5*(x1 + x2);
  double Dx = 0.5*(x2[0]-x1[0]);
  double Dy = 0.5*(x2[1]-x1[1]);
  double Dz = 0.5*(x2[2]-x1[2]);
  double D = max(Dx, max(Dy, Dz));
  x1 = xm - 2*vec3_t(D,D,D);
  x2 = xm + 2*vec3_t(D,D,D);
  m_OTGrid.setBounds(x1, x2);
  m_Length = (x2-x1).abs();
  m_OTGrid.setSmoothTransitionOff();
  m_OTGrid.setMaxCells(m_MaxOTCells);
  //EG_ERR_RETURN("stopped");
}

void SurfaceProjection::setBackgroundGrid_refineFromNodes()
{
  int num_refine;
  do {
    num_refine = 0;
    m_OTGrid.resetRefineMarks();
    foreach (vtkIdType id_node, m_Nodes) {
      vec3_t x;
      m_BGrid->GetPoints()->GetPoint(id_node, x.data());
      int i_otcell = m_OTGrid.findCell(x);
      double Dx = m_OTGrid.getDx(i_otcell);
      double Dy = m_OTGrid.getDy(i_otcell);
      double Dz = m_OTGrid.getDz(i_otcell);
      double D = max(Dx, max(Dy, Dz));
      if (D > m_EdgeLength[id_node]) {
        if (D > 0.5*m_MinOTLength) {
          m_OTGrid.markToRefine(i_otcell);
        }
      }
    }
    num_refine = m_OTGrid.refineAll();
    //cout << "refine from nodes: " << m_OTGrid.getNumCells() << "cells" << endl;
  } while (num_refine > 0);
}

void SurfaceProjection::setBackgroundGrid_refineFromEdges()
{
  int num_refine;
  do {
    num_refine = 0;
    m_OTGrid.resetRefineMarks();
    for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {
      vec3_t x1;
      m_BGrid->GetPoints()->GetPoint(m_Nodes[i_nodes], x1.data());
      for (int i_neigh = 0; i_neigh < m_N2N[i_nodes].size(); ++i_neigh) {
        if (i_nodes < i_neigh) {
          vec3_t x2;
          m_BGrid->GetPoints()->GetPoint(m_Nodes[i_neigh], x2.data());
          for (int i_cells = 0; i_cells < m_OTGrid.getNumCells(); ++i_cells) {
            if (!m_OTGrid.hasChildren(i_cells)) {
              double Dx = m_OTGrid.getDx(i_cells);
              double Dy = m_OTGrid.getDy(i_cells);
              double Dz = m_OTGrid.getDz(i_cells);
              double D = max(Dx, max(Dy, Dz));
              for (int i_faces = 0; i_faces < 6; ++i_faces) {
                double k;
                if (m_OTGrid.intersectsFace(i_cells, i_faces, x1, x2, k)) {
                  double L = min(m_EdgeLength[i_nodes], m_EdgeLength[i_neigh]); //(1-k)*m_EdgeLength[i_nodes] + k*m_EdgeLength[i_neigh];
                  if (D > L) {
                    if (D > 0.5*m_MinOTLength) {
                      m_OTGrid.markToRefine(i_cells);
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    num_refine = m_OTGrid.refineAll();
    //cout << "refine from edges: " << m_OTGrid.getNumCells() << "cells" << endl;
  } while (num_refine > 0);
}

void SurfaceProjection::setBackgroundGrid_refineFromFaces()
{
  int num_refine;
  do {
    num_refine = 0;
    m_OTGrid.resetRefineMarks();
    for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
      vtkIdType Npts, *pts;
      m_BGrid->GetCellPoints(id_cell, Npts, pts);
      if (Npts == 3) {
        vec3_t a, b, c;
        m_BGrid->GetPoints()->GetPoint(pts[0], a.data());
        m_BGrid->GetPoints()->GetPoint(pts[1], b.data());
        m_BGrid->GetPoints()->GetPoint(pts[2], c.data());
        double La = m_EdgeLength[pts[0]];
        double Lb = m_EdgeLength[pts[1]];
        double Lc = m_EdgeLength[pts[2]];
        for (int i_cells = 0; i_cells < m_OTGrid.getNumCells(); ++i_cells) {
          if (!m_OTGrid.hasChildren(i_cells) && !m_OTGrid.markedForRefine(i_cells)) {
            QVector<SortedPair<int> > edges;
            m_OTGrid.getEdges(i_cells, edges);
            foreach (SortedPair<int> edge, edges) {
              vec3_t xi;
              vec3_t ri;
              vec3_t x1 = m_OTGrid.getNodePosition(edge.v1);
              vec3_t x2 = m_OTGrid.getNodePosition(edge.v2);
              if (GeometryTools::intersectEdgeAndTriangle(a, b, c, x1, x2, xi, ri)) {
                double L = min(La, min(Lb, Lc)); //La + ri[0]*(Lb-La) + ri[1]*(Lc-La);
                double Dx = m_OTGrid.getDx(i_cells);
                double Dy = m_OTGrid.getDy(i_cells);
                double Dz = m_OTGrid.getDz(i_cells);
                double D = max(Dx, max(Dy, Dz));
                if (D > L) {
                  if (D > 0.5*m_MinOTLength) {
                    m_OTGrid.markToRefine(i_cells);
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
    num_refine = m_OTGrid.refineAll();
    //cout << "refine from faces: " << m_OTGrid.getNumCells() << "cells" << endl;
  } while (num_refine > 0);
}

void SurfaceProjection::setBackgroundGrid_computeLevelSet()
{
  // initialise G
  m_G.fill(0, m_OTGrid.getNumNodes());

  for (int i_nodes = 0; i_nodes < m_OTGrid.getNumNodes(); ++i_nodes) {
    double weight = 0;
    vec3_t xp = m_OTGrid.getNodePosition(i_nodes);
    foreach (Triangle T, m_Triangles) {
      vec3_t xi(1e99,1e99,1e99);
      vec3_t ri;
      double scal = (xp - T.m_a)*T.m_g3;
      vec3_t x1, x2;
      if (scal > 0) {
        x1 = xp + T.m_g3;
        x2 = xp - scal*T.m_g3 - T.m_g3;
      } else {
        x1 = xp - T.m_g3;
        x2 = xp - scal*T.m_g3 + T.m_g3;
      }
      double d = 1e99;

      bool intersects_face = GeometryTools::intersectEdgeAndTriangle(T.m_a, T.m_b, T.m_c, x1, x2, xi, ri);
      if (!intersects_face) {
        double kab = GeometryTools::intersection(T.m_a, T.m_b - T.m_a, xp, T.m_b - T.m_a);
        double kac = GeometryTools::intersection(T.m_a, T.m_c - T.m_a, xp, T.m_c - T.m_a);
        double kbc = GeometryTools::intersection(T.m_b, T.m_c - T.m_b, xp, T.m_c - T.m_b);
        double dab = (T.m_a + kab*(T.m_b-T.m_a) - xp).abs();
        double dac = (T.m_a + kac*(T.m_c-T.m_a) - xp).abs();
        double dbc = (T.m_b + kbc*(T.m_c-T.m_b) - xp).abs();
        bool set = false;
        if ((kab >= 0) && (kab <= 1)) {
          if (dab < d) {
            xi = T.m_a + kab*(T.m_b-T.m_a);
            d = dab;
            set = true;
          }
        }
        if ((kac >= 0) && (kac <= 1)) {
          if (dac < d) {
            xi = T.m_a + kac*(T.m_c-T.m_a);
            d = dac;
            set = true;
          }
        }
        if ((kbc >= 0) && (kbc <= 1)) {
          if (dbc < d) {
            xi = T.m_b + kbc*(T.m_c-T.m_b);
            d = dbc;
            set = true;
          }
        }
        double da = (T.m_a - xp).abs();
        double db = (T.m_b - xp).abs();
        double dc = (T.m_c - xp).abs();
        if (da < d) {
          xi = T.m_a;
          d = da;
          set = true;
        }
        if (db < d) {
          xi = T.m_b;
          d = db;
        }
        if (dc < d) {
          xi = T.m_c;
          d = dc;
          set = true;
        }
        if (!set) {
          EG_BUG;
        }
      }
      if (xi[0] > 1e98) {
        EG_BUG;
      }      
      double L = m_Length;
      vec3_t dx = xp - xi;
      double g = dx*T.m_g3;
      if (intersects_face) {
        d = fabs(g);
      }
      double w = 1;
      w *= m_DistWeight*pow(L/max(1e-6*L, d), m_DistExp);
      if (dx.abs() < 1e-6*L) {
        w *= m_DirWeight;
      } else {
        dx.normalise();
        w *= m_DirWeight*pow(fabs(dx*T.m_g3), m_DirExp);
      }
      w += m_WeightOffset;

      m_G[i_nodes] += w*g;
      weight += w;
    }
    m_G[i_nodes] /= weight;
  }

  // analytical test for sphere
  return;
  static bool first = true;
  if (first) {
    first = false;
    for (int i_nodes = 0; i_nodes < m_OTGrid.getNumNodes(); ++i_nodes) {
      vec3_t x = m_OTGrid.getNodePosition(i_nodes);
      double r = (x - vec3_t(-5,10,-5)).abs();
      m_G[i_nodes] = 1-r;
    }
  }
}

vec3_t SurfaceProjection::calcGradG(vec3_t x)
{
  int cell = m_OTGrid.findCell(x);
  vec3_t DG(0,0,0);
  DG[0] += m_G[m_OTGrid.getNode(cell, 1)] - m_G[m_OTGrid.getNode(cell, 0)];
  DG[0] += m_G[m_OTGrid.getNode(cell, 3)] - m_G[m_OTGrid.getNode(cell, 2)];
  DG[0] += m_G[m_OTGrid.getNode(cell, 5)] - m_G[m_OTGrid.getNode(cell, 4)];
  DG[0] += m_G[m_OTGrid.getNode(cell, 7)] - m_G[m_OTGrid.getNode(cell, 6)];
  DG[0] /= m_OTGrid.getDx(cell);
  DG[1] += m_G[m_OTGrid.getNode(cell, 2)] - m_G[m_OTGrid.getNode(cell, 0)];
  DG[1] += m_G[m_OTGrid.getNode(cell, 3)] - m_G[m_OTGrid.getNode(cell, 1)];
  DG[1] += m_G[m_OTGrid.getNode(cell, 6)] - m_G[m_OTGrid.getNode(cell, 4)];
  DG[1] += m_G[m_OTGrid.getNode(cell, 7)] - m_G[m_OTGrid.getNode(cell, 5)];
  DG[1] /= m_OTGrid.getDy(cell);
  DG[2] += m_G[m_OTGrid.getNode(cell, 4)] - m_G[m_OTGrid.getNode(cell, 0)];
  DG[2] += m_G[m_OTGrid.getNode(cell, 5)] - m_G[m_OTGrid.getNode(cell, 1)];
  DG[2] += m_G[m_OTGrid.getNode(cell, 6)] - m_G[m_OTGrid.getNode(cell, 2)];
  DG[2] += m_G[m_OTGrid.getNode(cell, 7)] - m_G[m_OTGrid.getNode(cell, 3)];
  DG[2] /= m_OTGrid.getDz(cell);
  DG *= 0.25;
  return DG;
}

double SurfaceProjection::calcG(vec3_t x)
{
  int cell = m_OTGrid.findCell(x);
  vec3_t r = x - m_OTGrid.getNodePosition(m_OTGrid.getNode(cell, 0));
  double kx = r[0]/m_OTGrid.getDx(cell);
  double ky = r[1]/m_OTGrid.getDy(cell);
  double kz = r[2]/m_OTGrid.getDz(cell);
  double g_01 = (1-kx)*m_G[m_OTGrid.getNode(cell, 0)] + kx*m_G[m_OTGrid.getNode(cell, 1)];
  double g_23 = (1-kx)*m_G[m_OTGrid.getNode(cell, 2)] + kx*m_G[m_OTGrid.getNode(cell, 3)];
  double g_45 = (1-kx)*m_G[m_OTGrid.getNode(cell, 4)] + kx*m_G[m_OTGrid.getNode(cell, 5)];
  double g_67 = (1-kx)*m_G[m_OTGrid.getNode(cell, 6)] + kx*m_G[m_OTGrid.getNode(cell, 7)];
  double g_01_23 = (1-ky)*g_01 + ky*g_23;
  double g_45_67 = (1-ky)*g_45 + ky*g_67;
  return (1-kz)*g_01_23 + kz*g_45_67;
}

void SurfaceProjection::writeOctree(QString file_name)
{
  EG_VTKSP(vtkUnstructuredGrid, otg);
  m_OTGrid.toVtkGrid(otg);
  EG_VTKSP(vtkDoubleArray, g);
  g->SetName("g");
  g->SetNumberOfValues(otg->GetNumberOfPoints());
  otg->GetPointData()->AddArray(g);
  for (int i = 0; i < otg->GetNumberOfPoints(); ++i) {
    g->SetValue(i, m_G[i]);
  }
  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
  vtu->SetFileName(qPrintable(GuiMainWindow::pointer()->getCwd() + "/" + file_name + ".vtu"));
  vtu->SetDataModeToBinary();
  vtu->SetInput(otg);
  vtu->Write();
  //writeGrid(m_BGrid, "m_BGrid");
}

vec3_t SurfaceProjection::projectWithLevelSet(vec3_t x)
{
  int count = 0;
  double g0 = calcG(x);
  double g = g0;
  if (g0 != 0) {
    while ((fabs(g/g0) > m_ConvLimit) && (count < m_MaxIter)) {
      ++count;
      vec3_t dx = calcGradG(x);
      dx.normalise();
      if (g != 0) {
        x -= m_Relax*g*dx;
      }
      g = calcG(x);
    }
  }
  return x;
}

void SurfaceProjection::writeGridWithNormals(QString filename)
{
  //qDebug()<<"void SurfaceProjection::writeGridWithNormals() called";
  
  vtkDoubleArray *vectors_normals = vtkDoubleArray::New();
  vectors_normals->SetName("normals");
  vectors_normals->SetNumberOfComponents(3);
  vectors_normals->SetNumberOfTuples(m_BGrid->GetNumberOfPoints());
  
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    vec3_t N = m_NodeNormals[id_node];
    double n[3];
    n[0]=N[0];
    n[1]=N[1];
    n[2]=N[2];
    vectors_normals->InsertTuple(id_node,n);
  }
  
  m_BGrid->GetPointData()->SetVectors(vectors_normals);
  
  vectors_normals->Delete();
  
/*  EG_VTKSP(vtkIntArray, var);
  var->SetName("");
  var->SetNumberOfValues(Ncells);
  grid->GetCellData()->AddArray(var);
  var->SetValue(i,0);
  
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    
  }*/
  
  saveGrid(m_BGrid, filename+"_BGrid_WithNormals");
}

void debug_output( QVector < QPair<vec3_t,vec3_t> > points,   QVector < QPair<vec3_t,vec3_t> > lines )
{
  QFile file("debug.vtk");
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    return;
  
  QTextStream out(&file);
  out<<"# vtk DataFile Version 2.0"<<endl;
  out<<"Unstructured Grid Example"<<endl;
  out<<"ASCII"<<endl;
  out<<""<<endl;
  out<<"DATASET UNSTRUCTURED_GRID"<<endl;
  out<<"POINTS "<<points.size()+2*lines.size()<<" double"<<endl;
  for(int i=0;i<points.size();i++) {
    vec3_t P = points[i].first;
    out<<P[0]<<" "<<P[1]<<" "<<P[2]<<endl;
  }
  for(int i=0;i<lines.size();i++) {
    vec3_t P1 = lines[i].first;
    vec3_t P2 = lines[i].second;
    out<<P1[0]<<" "<<P1[1]<<" "<<P1[2]<<endl;
    out<<P2[0]<<" "<<P2[1]<<" "<<P2[2]<<endl;
  }
  out<<""<<endl;
  out<<"CELLS "<<1+lines.size()<<" "<<4+3*lines.size()<<endl;
  out<<"3 0 1 2"<<endl;
  for(int i=0;i<lines.size();i++) {
    int P1 = points.size()+2*i;
    int P2 = P1+1;
    out<<2<<" "<<P1<<" "<<P2<<endl;
  }
  out<<""<<endl;
  out<<"CELL_TYPES "<<1+lines.size()<<endl;
  out<<"5"<<endl;
  for(int i=0;i<lines.size();i++) {
    out<<"3"<<endl;
  }
  out<<""<<endl;
  out<<"POINT_DATA "<<points.size()+2*lines.size()<<endl;
/*  out<<"SCALARS scalars double 1"<<endl;
  out<<"LOOKUP_TABLE default"<<endl;
  out<<"1.0"<<endl;
  out<<"2.0"<<endl;
  out<<"3.0"<<endl;*/
  out<<"VECTORS normals float"<<endl;
  for(int i=0;i<points.size();i++) {
    vec3_t N = points[i].second;
    out<<N[0]<<" "<<N[1]<<" "<<N[2]<<endl;
  }
  for(int i=0;i<lines.size();i++) {
    out<<0<<" "<<0<<" "<<0<<endl;
    out<<0<<" "<<0<<" "<<0<<endl;
  }
}

double interpolate(vec2_t M, vec2_t A, vec2_t nA, vec2_t B, vec2_t nB)
{
  double x0 = A[0];
  double x1 = B[0];
  double x = M[0];
  
  double f_x0 = A[1];
  double f_x1 = B[1];
  double fp_x0 = -nA[0]/nA[1];
  double fp_x1 = -nB[0]/nB[1];
  
  mat4_t matrix;
  matrix[0][0]=3*pow(x0,2); matrix[0][1]=2*x0;      matrix[0][2]=1;  matrix[0][3]=0;
  matrix[1][0]=3*pow(x1,2); matrix[1][1]=2*x1;      matrix[1][2]=1;  matrix[1][3]=0;
  matrix[2][0]=pow(x0,3);   matrix[2][1]=pow(x0,2); matrix[2][2]=x0; matrix[2][3]=1;
  matrix[3][0]=pow(x1,3);   matrix[3][1]=pow(x1,2); matrix[3][2]=x1; matrix[3][3]=1;
  
  mat4_t matrix_inv = matrix.inverse();
  vec4_t Y = vec4_t(fp_x0,fp_x1,f_x0,f_x1);
  vec4_t coeffs = matrix_inv * Y;
  
  // f(x)=a*x^3 + b*x^2 + c*x + d
  double ret = coeffs[0]*pow(x,3) + coeffs[1]*pow(x,2) + coeffs[2]*x + coeffs[3];
  return ret;
}

vec2_t interpolate_bezier(double t, vec2_t P0, vec2_t P1, vec2_t P2, vec2_t P3)
{
  // B(t)=(1-t^3)*P0 + 3*(1-t)^2*t*P1 + 3*(1-t)*t^2*P2 + t^3*P3;
  return (1-pow(t,3))*P0 + 3*pow((1-t),2)*t*P1 + 3*(1-t)*pow(t,2)*P2 + pow(t,3)*P3;
}

vec3_t interpolate_bezier(double t, vec3_t P0, vec3_t P1, vec3_t P2, vec3_t P3)
{
  // B(t)=(1-t^3)*P0 + 3*(1-t)^2*t*P1 + 3*(1-t)*t^2*P2 + t^3*P3;
  return (1-pow(t,3))*P0 + 3*pow((1-t),2)*t*P1 + 3*(1-t)*pow(t,2)*P2 + pow(t,3)*P3;
}

vec2_t interpolate_bezier_normals(double t, vec2_t P0, vec2_t N0, vec2_t P3, vec2_t N3)
{
  vec2_t P1(N0[1],-N0[0]);
  vec2_t P2(-N3[1],N3[0]);
  return interpolate_bezier(t,P0,P1,P2,P3);
}

vec3_t interpolate_bezier_normals(double t, vec3_t P0, vec3_t N0, vec3_t P3, vec3_t N3, vec3_t u, vec3_t v)
{
  double N0x=N0*u/u.abs2();
  double N0y=N0*v/v.abs2();
  double N3x=N3*u/u.abs2();
  double N3y=N3*v/v.abs2();
  
  vec3_t P0P1 = (N0y)*u + (-N0x)*v;
  vec3_t P3P2 = (-N3y)*u + (N3x)*v;
  
  P0P1.normalise();
  P3P2.normalise();
  
  vec3_t P1 = P0 + P0P1;
  vec3_t P2 = P3 + P3P2;
  return interpolate_bezier(t,P0,P1,P2,P3);
}

// coordinate systems:
// 3D:
// global: X,Y,Z -> g_
// local: g1,g2,g3 -> l_
// 2D:
// triangle: g1,g2 -> t_
// middle planes
// plane 1: AI1, g3 -> pm1_
// plane 2: BI2, g3 -> pm2_
// plane 3: CI3, g3 -> pm3_
// orthogonal edge planes
// plane 1: BC, g3 -> poe1_
// plane 2: CA, g3 -> poe2_
// plane 3: AB, g3 -> poe3_
// non-orthogonal edge planes
// plane 1: BC, nI3 -> pnoe1_
// plane 2: CA, nI1 -> pnoe2_
// plane 3: AB, nI2 -> pnoe3_
vec3_t SurfaceProjection::correctCurvature1(int i_tri, vec3_t g_M)
{
//   qDebug()<<"vec3_t SurfaceProjection::correctCurvature1 called";
  
  Triangle T = m_Triangles[i_tri];
  vec3_t g_A = T.m_a;
  vec3_t g_B = T.m_b;
  vec3_t g_C = T.m_c;
  
  vec3_t l_A(0,0,0);
  vec3_t l_B(1,0,0);
  vec3_t l_C(0,1,0);
  vec3_t l_M = T.global3DToLocal3D(g_M);
  
  vec2_t t_A(0,0);
  vec2_t t_B(1,0);
  vec2_t t_C(0,1);
  vec2_t t_M(l_M[0],l_M[1]);
  
  vec3_t g_nA = m_NodeNormals[T.m_id_a];
  vec3_t g_nB = m_NodeNormals[T.m_id_b];
  vec3_t g_nC = m_NodeNormals[T.m_id_c];
  
  vec2_t pm1_A(0,0);
  vec2_t pm1_I1(1,0);
  
  vec2_t pm2_B(0,0);
  vec2_t pm2_I2(1,0);
  
  vec2_t pm3_C(0,0);
  vec2_t pm3_I3(1,0);
  
  vec3_t g_g1 = T.m_g1;
  vec3_t g_g2 = T.m_g2;
  vec3_t g_g3 = T.m_g3;
  
  vec3_t l_g1(1,0,0);
  vec3_t l_g2(0,1,0);
  vec3_t l_g3(0,0,1);
  
  // calculating intersections
  double k1,k2;
  if(!intersection (k1, k2, t_A, t_M-t_A, t_B, t_C-t_B)) return(g_M);
  vec2_t t_I1 = t_A+k1*(t_M-t_A);
  vec3_t g_nI1 = (1-k2)*g_nB + k2*g_nC;
  vec2_t pm1_M(1.0/k1,0);
  vec2_t poe1_I1(k2,0);
  vec2_t pnoe1_I1(k2,0);
  
  if(!intersection (k1, k2, t_B, t_M-t_B, t_C, t_A-t_C)) return(g_M);
  vec2_t t_I2 = t_B+k1*(t_M-t_B);
  vec3_t g_nI2 = (1-k2)*g_nC + k2*g_nA;
  vec2_t pm2_M(1.0/k1,0);
  vec2_t poe2_I2(k2,0);
  vec2_t pnoe2_I2(k2,0);
  
  if(!intersection (k1, k2, t_C, t_M-t_C, t_A, t_B-t_A)) return(g_M);
  vec2_t t_I3 = t_C+k1*(t_M-t_C);
  vec3_t g_nI3 = (1-k2)*g_nA + k2*g_nB;
  vec2_t pm3_M(1.0/k1,0);
  vec2_t poe3_I3(k2,0);
  vec2_t pnoe3_I3(k2,0);
  
  // after knowing intersections
  vec3_t l_I1(t_I1[0],t_I1[1],0);
  vec3_t l_I2(t_I2[0],t_I2[1],0);
  vec3_t l_I3(t_I3[0],t_I3[1],0);
  
  vec3_t g_I1 = g_A+T.m_G*l_I1;
  vec3_t g_I2 = g_A+T.m_G*l_I2;
  vec3_t g_I3 = g_A+T.m_G*l_I3;
  
  vec3_t l_nI1 = T.m_GI*g_nI1;
  vec3_t l_nI2 = T.m_GI*g_nI2;
  vec3_t l_nI3 = T.m_GI*g_nI3;
  
  vec3_t l_nA = T.m_GI*g_nA;
  vec3_t l_nB = T.m_GI*g_nB;
  vec3_t l_nC = T.m_GI*g_nC;
  
  vec3_t l_AI1 = l_I1 - l_A;
  vec3_t l_BI2 = l_I2 - l_B;
  vec3_t l_CI3 = l_I3 - l_C;
  
  vec3_t g_AI1 = g_I1 - g_A;
  vec3_t g_BI2 = g_I2 - g_B;
  vec3_t g_CI3 = g_I3 - g_C;
  
  vec2_t pm1_nA = projectVectorOnPlane(l_nA,l_AI1,l_g3);
  vec2_t pm1_nI1 = projectVectorOnPlane(l_nI1,l_AI1,l_g3);
  
  vec2_t pm2_nB = projectVectorOnPlane(l_nB,l_BI2,l_g3);
  vec2_t pm2_nI2 = projectVectorOnPlane(l_nI2,l_BI2,l_g3);
  
  vec2_t pm3_nC = projectVectorOnPlane(l_nC,l_CI3,l_g3);
  vec2_t pm3_nI3 = projectVectorOnPlane(l_nI3,l_CI3,l_g3);
  
  // plane 1: BC, nI1 -> pnoe1_
  // plane 2: CA, nI2 -> pnoe2_
  // plane 3: AB, nI3 -> pnoe3_
  
  vec3_t l_AB = l_B - l_A;
  vec3_t l_BC = l_C - l_B;
  vec3_t l_CA = l_A - l_C;
  
  ///////////////
  vec2_t pnoe1_B(0,0);
  vec2_t pnoe1_C(1,0);
  vec2_t pnoe1_nB = projectVectorOnPlane(l_nB,l_BC,l_nI1);
  vec2_t pnoe1_nC = projectVectorOnPlane(l_nC,l_BC,l_nI1);
  vec2_t pnoe1_nI1(0,1);
  vec2_t pnoe1_tB = turnRight(pnoe1_nB);
  vec2_t pnoe1_tC = turnRight(pnoe1_nC);
  ///////////////
  vec2_t pnoe2_C(0,0);
  vec2_t pnoe2_A(1,0);
  vec2_t pnoe2_nC = projectVectorOnPlane(l_nC,l_CA,l_nI2);
  vec2_t pnoe2_nA = projectVectorOnPlane(l_nA,l_CA,l_nI2);
  vec2_t pnoe2_nI2(0,1);
  vec2_t pnoe2_tC = turnRight(pnoe2_nC);
  vec2_t pnoe2_tA = turnRight(pnoe2_nA);
  ///////////////
  vec2_t pnoe3_A(0,0);
  vec2_t pnoe3_B(1,0);
  vec2_t pnoe3_nA = projectVectorOnPlane(l_nA,l_AB,l_nI3);
  vec2_t pnoe3_nB = projectVectorOnPlane(l_nB,l_AB,l_nI3);
  vec2_t pnoe3_nI3(0,1);
  vec2_t pnoe3_tA = turnRight(pnoe3_nA);
  vec2_t pnoe3_tB = turnRight(pnoe3_nB);
  ///////////////
  vec3_t g_J1;
  vec3_t g_J2;
  vec3_t g_J3;
  
  
  // interpolation attempts
  double z1 = interpolate(pm1_M, pm1_A, pm1_nA, pm1_I1, pm1_nI1);
  double z2 = interpolate(pm2_M, pm2_B, pm2_nB, pm2_I2, pm2_nI2);
  double z3 = interpolate(pm3_M, pm3_C, pm3_nC, pm3_I3, pm3_nI3);
  double z = z1;//(z1+z2+z3)/3.0;
  
  vec3_t l_X = l_M + z*l_g3;
  vec3_t g_X = g_A+T.m_G*l_X;
  
  // returning value
  return g_X;
}// end of correctCurvature

vec3_t SurfaceProjection::cylinder(vec3_t center, double radius, vec3_t g_M)
{
  vec3_t g_P;
  g_P[2]=g_M[2];
  double L = sqrt(g_M[0]*g_M[0]+g_M[1]*g_M[1]);
  g_P[0]=radius * g_M[0]/L;
  g_P[1]=radius * g_M[1]/L;
  return g_P;
}

vec3_t SurfaceProjection::cylinder(vec3_t center, double radius, int i_tri, vec3_t r)
{
  Triangle T = m_Triangles[i_tri];
  vec3_t g_A = T.m_a;
  vec3_t g_M = g_A+T.m_G*r;
  return cylinder(center, radius, g_M);
}

vec3_t SurfaceProjection::projectWithGeometry_original(vec3_t xp, vtkIdType id_node)
{
  vec3_t x_proj(1e99,1e99,1e99), r_proj;
  bool on_triangle = false;
  bool need_full_search = false;
  if (id_node >= m_ProjTriangles.size()) {
    int old_size = m_ProjTriangles.size();
    m_ProjTriangles.resize(m_FGrid->GetNumberOfPoints());
    for (int i = old_size; i < m_ProjTriangles.size(); ++i) {
      m_ProjTriangles[i] = -1;
    }
    need_full_search = true;
  } else {
    int i_triangles = m_ProjTriangles[id_node];
    if (i_triangles < 0) {
      need_full_search = true;
    } else {
      if (i_triangles >= m_Triangles.size()) {
        EG_BUG;
      }
      Triangle T = m_Triangles[i_triangles];
      vec3_t xi, ri;
      double d;
      int side;
      bool intersects = T.projectOnTriangle(xp, xi, ri, d, side, true);
      if (!intersects || (d > 0.1*T.m_smallest_length)) {
        need_full_search = true;
      } else {
        x_proj = xi;
        r_proj = ri;
        on_triangle = intersects;
        ++m_NumDirect;
      }
    }
  }
  if (need_full_search) {
    ++m_NumFull;
    double d_min = 1e99;
    for (int i_triangles = 0; i_triangles < m_Triangles.size(); ++i_triangles) {
      Triangle T = m_Triangles[i_triangles];
      double d;
      vec3_t xi, ri;
      int side;
      bool intersects = T.projectOnTriangle(xp, xi, ri, d, side, true);
      if (d < d_min) {
        x_proj = xi;
        r_proj = ri;
        d_min = d;
        m_ProjTriangles[id_node] = i_triangles;
        on_triangle = intersects;
      }
    }
    if (x_proj[0] > 1e98) {
      EG_BUG;
    }
  }
  if (on_triangle) {
    //x_proj = correctCurvature(m_ProjTriangles[id_node], r_proj);
  }
  return x_proj;
}

vec3_t SurfaceProjection::projectWithGeometry(vec3_t xp, vtkIdType id_node)
{
  checkVector(xp);
  
/*  if(m_ExactMode==1) return ellipsoid(xp);
  if(m_ExactMode==2) return ellipse(xp);
  if(m_ExactMode==3) return rectangle(xp);
  if(m_ExactMode==4) return cuboid(xp);
  if(m_ExactMode==5) return cylinder(xp);*/
  
//   qDebug()<<"=== m_correctCurvature="<<m_correctCurvature<<" ===";
  
//   qWarning()<<"@@@@@@@@@@@@ xp="<<xp[0]<<xp[1]<<xp[2]<<endl;
  vec3_t x_proj(1e99,1e99,1e99), r_proj(0,0,0);
  bool x_proj_set =false;
  
  // initilizing booleans
  bool on_triangle = false;
  bool need_full_search = false;
  
  if (id_node >= m_ProjTriangles.size()) { //if there is no known triangle on which to project
    int old_size = m_ProjTriangles.size();
    m_ProjTriangles.resize(m_FGrid->GetNumberOfPoints());
    for (int i = old_size; i < m_ProjTriangles.size(); ++i) {
      m_ProjTriangles[i] = -1;
    }
    need_full_search = true;
  } else { //if there is a known triangle on which to project
    int i_triangles = m_ProjTriangles[id_node];
    if (i_triangles < 0) {
      need_full_search = true;
    } else {
      if (i_triangles >= m_Triangles.size()) {
        EG_BUG;
      }
      Triangle T = m_Triangles[i_triangles];
      vec3_t xi, ri;
      double d;
      int side;
      bool intersects = T.projectOnTriangle(xp, xi, ri, d, side, true);
      if (!intersects || (d > 0.1*T.m_smallest_length)) {
        need_full_search = true;
      } else {
        x_proj = xi; x_proj_set = true;
        if ( x_proj[0]>1e98 ) { // should never happen
          EG_BUG;
        }
        r_proj = ri;
        on_triangle = intersects;
        ++m_NumDirect;
      }
    }
  }
  if (true || need_full_search) {
//     qDebug()<<"starting full search";
    ++m_NumFull;
    double d_min = 1e99;
    bool first = true;
    for (int i_triangles = 0; i_triangles < m_Triangles.size(); ++i_triangles) {
      Triangle T = m_Triangles[i_triangles];
      double d;
      vec3_t xi, ri;
      int side;
      bool intersects = T.projectOnTriangle(xp, xi, ri, d, side, true);
//       if(d>9000) qWarning()<<"d="<<d;
      if(d>=1e99) EG_BUG;
      if (/*first ||*/ d < d_min) {
        x_proj = xi; x_proj_set = true;
        if ( x_proj[0]>1e98 ) { // should never happen
          EG_BUG;
        }
        r_proj = ri;
        d_min = d;
        m_ProjTriangles[id_node] = i_triangles;
        on_triangle = intersects;
        first = false;
      }
//       qDebug()<<"full search done";
    }
    if ( !x_proj_set ) { // should never happen
      checkVector(xp);
      qWarning()<<"No projection found for point xp="<<xp[0]<<xp[1]<<xp[2]<<endl;
      writeGrid(GuiMainWindow::pointer()->getGrid(),"griddump");
      EG_BUG;
    }
  }
//    if(on_triangle) {
//      if(m_correctCurvature) x_proj = correctCurvature(m_ProjTriangles[id_node], r_proj);
     if(m_correctCurvature) x_proj = correctCurvature2(m_ProjTriangles[id_node], xp);
//    }
  if(!on_triangle) {
//     qDebug()<<"WARNING: Not on triangle! id_node="<<id_node;
  }
  
//   writeGrid(m_BGrid,"m_BGrid");
  
  return x_proj;
}

vec3_t SurfaceProjection::project(vec3_t x, vtkIdType id_node)
{
  checkVector(x);
  
  if (m_UseLevelSet) {
    x = projectWithLevelSet(x);
  } else {
    if (id_node < 0) {
      EG_BUG;
    }
    x = projectWithGeometry_original(x, id_node);
  }
//   writeGrid(m_FGrid,"m_FGrid");
//   writeGrid(m_BGrid,"m_BGrid");
  return x;
}

int SurfaceProjection::getControlPoints_orthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110, double Lmax)
{
  vec3_t A = T.m_a;
  vec3_t B = T.m_b;
  vec3_t C = T.m_c;
  vec3_t nA = m_NodeNormals[T.m_id_a];
  vec3_t nB = m_NodeNormals[T.m_id_b];
  vec3_t nC = m_NodeNormals[T.m_id_c];
  
  //cout<<"nA="<<nA<<endl;
  //cout<<"nB="<<nB<<endl;
  //cout<<"nC="<<nC<<endl;

//   cout<<"A="<<A<<" B="<<B<<" C="<<C<<endl;
  //cout<<"-->BC"<<endl;
  
  if(T.m_Valid) {
    X_011 = intersectionOnPlane(T.m_g3, B, nB, C, nC);
    X_101 = intersectionOnPlane(T.m_g3, C, nC, A, nA);
    X_110 = intersectionOnPlane(T.m_g3, A, nA, B, nB);
  }
  else {
    X_011 = 0.5*(B+C);
    X_101 = 0.5*(C+A);
    X_110 = 0.5*(A+B);
  }
  
  if (!checkVector(X_011)) EG_BUG;
  if (!checkVector(X_101)) EG_BUG;
  if (!checkVector(X_110)) EG_BUG;
  
  limitControlPoints(T, X_011, X_101, X_110);
  
  if (!checkVector(X_011)) EG_BUG;
  if (!checkVector(X_101)) EG_BUG;
  if (!checkVector(X_110)) EG_BUG;
  
  return(0);
}

int SurfaceProjection::getControlPoints_nonorthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110, double Lmax)
{
  vec3_t A = T.m_a;
  vec3_t B = T.m_b;
  vec3_t C = T.m_c;
  vec3_t nA = m_NodeNormals[T.m_id_a];
  vec3_t nB = m_NodeNormals[T.m_id_b];
  vec3_t nC = m_NodeNormals[T.m_id_c];
  
//   cout<<"A="<<A<<" B="<<B<<" C="<<C<<endl;
  
  if ((0.5*(nB+nC)).abs2()==0) EG_BUG;
  if ((0.5*(nC+nA)).abs2()==0) EG_BUG;
  if ((0.5*(nA+nB)).abs2()==0) EG_BUG;
  
  if(T.m_Valid) {
    X_011 = intersectionOnPlane(0.5*(nB+nC), B, nB, C, nC);
    X_101 = intersectionOnPlane(0.5*(nC+nA), C, nC, A, nA);
    X_110 = intersectionOnPlane(0.5*(nA+nB), A, nA, B, nB);
  }
  else {
    X_011 = 0.5*(B+C);
    X_101 = 0.5*(C+A);
    X_110 = 0.5*(A+B);
  }
  
  /// \todo make sure nBC,nCA,nAB are not null vectors!!!
  vec3_t nBC = 0.5*(nB+nC);
  vec3_t nCA = 0.5*(nC+nA);
  vec3_t nAB = 0.5*(nA+nB);
  
  if (!checkVector(X_011)) {
    qWarning()<<X_011<<" = intersectionOnPlane("<<nBC<<", "<<B<<", "<<nB<<", "<<C<<", "<<nC<<")";
    EG_BUG;
  }
  if (!checkVector(X_101)) {
    qWarning()<<X_101<<" = intersectionOnPlane("<<nCA<<", "<<C<<", "<<nC<<", "<<A<<", "<<nA<<")";
    EG_BUG;
  }
  if (!checkVector(X_110)) {
    qWarning()<<X_110<<" = intersectionOnPlane("<<nAB<<", "<<A<<", "<<nA<<", "<<B<<", "<<nB<<")";
    EG_BUG;
  }
  
  limitControlPoints(T, X_011, X_101, X_110);
  
  if (!checkVector(X_011)) EG_BUG;
  if (!checkVector(X_101)) EG_BUG;
  if (!checkVector(X_110)) EG_BUG;
  
  return(0);
}

int SurfaceProjection::limitControlPoints(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110)
{
  vec3_t A=T.m_a;
  vec3_t B=T.m_b;
  vec3_t C=T.m_c;
  vec3_t nA = m_NodeNormals[T.m_id_a];
  vec3_t nB = m_NodeNormals[T.m_id_b];
  vec3_t nC = m_NodeNormals[T.m_id_c];
  
/*  vec3_t P_011 = projectPointOnEdge(X_011,B,C-B);
  vec3_t P_101 = projectPointOnEdge(X_101,C,A-C);
  vec3_t P_110 = projectPointOnEdge(X_110,A,B-A);*/
  
  vec3_t P_011 = 0.5*(B+C);
  vec3_t P_101 = 0.5*(A+C);
  vec3_t P_110 = 0.5*(A+B);
  
//   Lmax = 1.0*T.m_smallest_length;
  
  double Lmax_011 = (C - B).abs();
  double Lmax_101 = (A - C).abs();
  double Lmax_110 = (B - A).abs();
  
  double L_011 = (X_011-P_011).abs();
  double L_101 = (X_101-P_101).abs();
  double L_110 = (X_110-P_110).abs();
  
  checkVector(X_011);
  checkVector(P_011);
  checkVector(X_101);
  checkVector(P_101);
  checkVector(X_110);
  checkVector(P_110);
  
  if( L_011 > Lmax_011 ) {
//     qWarning()<<"WARNING: CONTROL POINT RESTRICTED: Lmax="<<Lmax;
//     qWarning()<<"X_011="<<X_011<<"P_011="<<P_011<<"L_011="<<L_011;
    X_011 = P_011 + Lmax_011/L_011 * (X_011-P_011);
  }
  if( L_101 > Lmax_101 ) {
//     qWarning()<<"WARNING: CONTROL POINT RESTRICTED: Lmax="<<Lmax;
//     qWarning()<<"X_101="<<X_101<<"P_101="<<P_101<<"L_101="<<L_101;
    X_101 = P_101 + Lmax_101/L_101 * (X_101-P_101);
  }
  if( L_110 > Lmax_110 ) {
//     qWarning()<<"WARNING: CONTROL POINT RESTRICTED: Lmax="<<Lmax;
//     qWarning()<<"X_110="<<X_110<<"P_110="<<P_110<<"L_110="<<L_110;
    X_110 = P_110 + Lmax_110/L_110 * (X_110-P_110);
  }
  return(0);
}

void SurfaceProjection::writeTriangleGrid(QString filename)
{
  int N_cells = m_Triangles.size();
  int N_points = 3*m_Triangles.size();
  
  qDebug()<<"N_cells="<<N_cells;
  qDebug()<<"N_points="<<N_points;
  
  
  EG_VTKSP(vtkUnstructuredGrid, triangle_grid);
  allocateGrid(triangle_grid , N_cells, N_points);
  
  vtkDoubleArray *vectors_hasneighbour0 = vtkDoubleArray::New();
  vectors_hasneighbour0->SetName("hasneighbour0");
  vectors_hasneighbour0->SetNumberOfComponents(3);
  vectors_hasneighbour0->SetNumberOfTuples(N_points);
  
  vtkDoubleArray *vectors_hasneighbour1 = vtkDoubleArray::New();
  vectors_hasneighbour1->SetName("hasneighbour1");
  vectors_hasneighbour1->SetNumberOfComponents(3);
  vectors_hasneighbour1->SetNumberOfTuples(N_points);
  
  vtkDoubleArray *vectors_hasneighbour2 = vtkDoubleArray::New();
  vectors_hasneighbour2->SetName("hasneighbour2");
  vectors_hasneighbour2->SetNumberOfComponents(3);
  vectors_hasneighbour2->SetNumberOfTuples(N_points);
  
  vtkDoubleArray *vectors_neighbours = vtkDoubleArray::New();
  vectors_neighbours->SetName("neighbours");
  vectors_neighbours->SetNumberOfComponents(3);
  vectors_neighbours->SetNumberOfTuples(N_cells);
  
  int node_count = 0;
  int cell_count = 0;
  for(int i=0; i<m_Triangles.size(); i++) {
    vtkIdType pts[3];
    triangle_grid->GetPoints()->SetPoint(node_count, m_Triangles[i].m_a.data()); pts[0]=node_count; node_count++;
    triangle_grid->GetPoints()->SetPoint(node_count, m_Triangles[i].m_b.data()); pts[1]=node_count; node_count++;
    triangle_grid->GetPoints()->SetPoint(node_count, m_Triangles[i].m_c.data()); pts[2]=node_count; node_count++;
  
    vec3_t v0,v1,v2;
    v0 = getEdgeNormal(m_Triangles[i].m_id_a,m_Triangles[i].m_id_b).normalise();
    v1 = getEdgeNormal(m_Triangles[i].m_id_b,m_Triangles[i].m_id_c).normalise();
    v2 = getEdgeNormal(m_Triangles[i].m_id_c,m_Triangles[i].m_id_a).normalise();
    if(!checkVector(v0) || !checkVector(v1) || !checkVector(v2)) {
      qWarning()<<"v0="<<v0;
      qWarning()<<"v1="<<v1;
      qWarning()<<"v2="<<v2;
      EG_BUG;
    };
    vectors_hasneighbour0->InsertTuple(pts[0],v0.data());
    vectors_hasneighbour1->InsertTuple(pts[1],v1.data());
    vectors_hasneighbour2->InsertTuple(pts[2],v2.data());
  
    vec3_t neighbours = m_Triangles[i].m_g3;
    if(!m_Triangles[i].m_has_neighbour[0]) neighbours+=v0;
    if(!m_Triangles[i].m_has_neighbour[1]) neighbours+=v1;
    if(!m_Triangles[i].m_has_neighbour[2]) neighbours+=v2;
    
    int cell = triangle_grid->InsertNextCell(VTK_TRIANGLE,3,pts);cell_count++;
    vectors_neighbours->InsertTuple(cell, neighbours.data());
    
  }
  
  triangle_grid->GetPointData()->AddArray(vectors_hasneighbour0);
  triangle_grid->GetPointData()->AddArray(vectors_hasneighbour1);
  triangle_grid->GetPointData()->AddArray(vectors_hasneighbour2);
  triangle_grid->GetCellData()->AddArray(vectors_neighbours);
  vectors_hasneighbour0->Delete();
  vectors_hasneighbour1->Delete();
  vectors_hasneighbour2->Delete();
  vectors_neighbours->Delete();
  
  saveGrid(triangle_grid, filename+"_TriangleGrid");
  
}

void SurfaceProjection::writeInterpolationGrid(QString filename)
{
  int N_cells = m_BGrid->GetNumberOfCells()+2*m_BGrid->GetNumberOfCells();
  int N_points = m_BGrid->GetNumberOfPoints()+6*m_BGrid->GetNumberOfCells();
  int N = 10;
  int N_cells_per_triangle = (N-1)*(N-1);
  int N_points_per_triangle = (N*N+N)/2;
  
  //qDebug()<<"N_cells="<<N_cells;
  //qDebug()<<"N_points="<<N_points;
  //qDebug()<<"N_cells_per_triangle="<<N_cells_per_triangle;
  //qDebug()<<"N_points_per_triangle="<<N_points_per_triangle;
  
  allocateGrid(m_InterpolationGrid , N_cells, N_points);
  makeCopyNoAlloc(m_BGrid, m_InterpolationGrid);
  
  MeshPartition new_grid_partition;
  EG_VTKSP(vtkUnstructuredGrid, bezier_first);
  bool first = true;
  
//   allocateGrid(m_BezierGrid, m_Triangles.size()*N_cells_per_triangle, m_Triangles.size()*N_points_per_triangle);

  vtkIdType node_count = m_BGrid->GetNumberOfPoints();
  int cell_count = m_BGrid->GetNumberOfCells();
  
  vtkIdType offset = 0;
  
  for (int i_triangles = 0; i_triangles < m_Triangles.size(); ++i_triangles) {
    Triangle T = m_Triangles[i_triangles];
    if ( i_triangles == 1 ) {
      //cout<<"+++++++++++++++++++++++++"<<endl;
    }
    vec3_t J1,K1;
    vec3_t J2,K2;
    vec3_t J3,K3;
    //qDebug()<<"=== ORTHOGONAL PLANES ===";
    getControlPoints_orthogonal(T,J1,J2,J3, 1e99);
    //qDebug()<<"=== NON-ORTHOGONAL PLANES ===";
    getControlPoints_nonorthogonal(T,K1,K2,K3, 1e99);
    
    vtkIdType idx_J1, idx_J2, idx_J3;
    m_InterpolationGrid->GetPoints()->SetPoint(node_count, J1.data()); idx_J1=node_count; node_count++;
    m_InterpolationGrid->GetPoints()->SetPoint(node_count, J2.data()); idx_J2=node_count; node_count++;
    m_InterpolationGrid->GetPoints()->SetPoint(node_count, J3.data()); idx_J3=node_count; node_count++;
    vtkIdType idx_K1, idx_K2, idx_K3;
    m_InterpolationGrid->GetPoints()->SetPoint(node_count, K1.data()); idx_K1=node_count; node_count++;
    m_InterpolationGrid->GetPoints()->SetPoint(node_count, K2.data()); idx_K2=node_count; node_count++;
    m_InterpolationGrid->GetPoints()->SetPoint(node_count, K3.data()); idx_K3=node_count; node_count++;
    
    if ( i_triangles == 1 ) {
      //cout<<"+++++++++++++++++++++++++"<<endl;
      //cout<<"A="<<T.m_a<<" B="<<T.m_b<<" C="<<T.m_c<<endl;
      //cout<<"J1="<<J1<<" K1="<<K1<<endl;
      //cout<<"J2="<<J2<<" K2="<<K2<<endl;
      //cout<<"J3="<<J3<<" K3="<<K3<<endl;
      //cout<<"+++++++++++++++++++++++++"<<endl;
    }
    
    // add the local grid
    BezierTriangle bezier_triangle(T.m_a, T.m_b, T.m_c, K1, K2, K3);
    if(first) {
      first = false;
      bezier_triangle.getBezierSurface(bezier_first, N);
      new_grid_partition.setGrid(bezier_first);
      new_grid_partition.setAllCells();
    }
    else {
      EG_VTKSP(vtkUnstructuredGrid, bezier);
      bezier_triangle.getBezierSurface(bezier, N);
      MeshPartition grid_partition(bezier, true);
      new_grid_partition.addPartition(grid_partition);
    }
    
    vtkIdType polyline_ortho[7];
    vtkIdType polyline_nonortho[7];
    
    polyline_ortho[0]=T.m_id_a;
    polyline_ortho[1]=idx_J3;
    polyline_ortho[2]=T.m_id_b;
    polyline_ortho[3]=idx_J1;
    polyline_ortho[4]=T.m_id_c;
    polyline_ortho[5]=idx_J2;
    polyline_ortho[6]=T.m_id_a;
    
    polyline_nonortho[0]=T.m_id_a;
    polyline_nonortho[1]=idx_K3;
    polyline_nonortho[2]=T.m_id_b;
    polyline_nonortho[3]=idx_K1;
    polyline_nonortho[4]=T.m_id_c;
    polyline_nonortho[5]=idx_K2;
    polyline_nonortho[6]=T.m_id_a;
    
    m_InterpolationGrid->InsertNextCell(4,7,polyline_ortho);cell_count++;
    m_InterpolationGrid->InsertNextCell(4,7,polyline_nonortho);cell_count++;
    
  }
  
  //qDebug()<<"node_count="<<node_count;
  //qDebug()<<"cell_count="<<cell_count;
  //qDebug()<<"offset="<<offset;

  saveGrid(m_InterpolationGrid, filename+"_InterpolationGrid");
  makeCopy(new_grid_partition.getGrid(), m_BezierGrid);
  saveGrid(m_BezierGrid, filename+"_BezierGrid");
//   saveGrid(new_grid_partition.getGrid(), filename+"_BezierGrid");
  this->writeGrid(m_BGrid,filename+"_BGrid");
  
}

// mat2_t SurfaceProjection::Jacobian_Matrix(double x, double y)
// {
//   mat2_t J;
//   return J;
// }

// vec2_t BezierProjectionFunction(double x, double y)
// {
//   vec2_t F;
//   quadraticBezierTriangle(double u, double v, double w, vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110);
//   return F;
// }

vec3_t SurfaceProjection::getEdgeNormal(vtkIdType id_node1, vtkIdType id_node2)
{
  l2l_t  n2n   = getPartN2N();
  l2l_t  n2c = getPartN2C();
  l2g_t cells = getPartCells();
  g2l_t _cells = getPartLocalCells();
  l2g_t nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  
  vec3_t x1,x2;
  m_Grid->GetPoints()->GetPoint(id_node1, x1.data());
  m_Grid->GetPoints()->GetPoint(id_node2, x2.data());
  
  QVector <vtkIdType> edge_cells;
  getEdgeCells(id_node1,id_node2,edge_cells);
  vtkIdType id_cell = edge_cells[0];
  int side = getSide(id_cell,m_BGrid,id_node1,id_node2);
  vtkIdType *pts, N_pts;
  m_BGrid->GetCellPoints(id_cell, N_pts, pts);
  vec3_t t;
  if(pts[side]==id_node1) {
    t = x2-x1;
  }
  else {
    t = x1-x2;
  }
  vec3_t n = m_Triangles[id_cell].m_g3;
  vec3_t Nedge = t.cross(n);
//   qDebug()<<"Nedge="<<Nedge;
  return Nedge;
}

vec3_t SurfaceProjection::ellipsoid(vec3_t M)
{
  vec3_t A = m_center;
  vec3_t AM = M-A;
  return( A + m_Rx.abs() * AM.normalise() );
}

vec3_t SurfaceProjection::ellipse(vec3_t M)
{
  vec3_t A = m_center;
  vec3_t AM = M-A;
  double x = AM*m_Rx;
  double y = AM*m_Ry;
  double z = AM*m_Rz;
//   qWarning()<<x<<y<<z;
  vec3_t P = A + x*m_Rx + y*m_Ry;
//   qWarning()<<P;
  vec3_t AP = P-A;
//   qWarning()<<AP;
//   qWarning()<<AP.normalise();
//   qWarning()<<m_Rx.abs();
//   qWarning()<<A + m_Rx.abs() * AP.normalise();
  if(AP.abs()<=m_Rx.abs()) return P;
  else return( A + m_Rx.abs() * AP.normalise() );
}

vec3_t SurfaceProjection::rectangle(vec3_t M)
{
  return vec3_t(0,0,0);
}

vec3_t SurfaceProjection::cuboid(vec3_t M)
{
  return vec3_t(0,0,0);
}

vec3_t SurfaceProjection::cylinder(vec3_t M)
{
  vec3_t A = m_center;
  vec3_t u = m_Rz;
  double k = (M-A)*u;
  vec3_t P = A + k*u;
  vec3_t PM = M-P;
  return P + m_Rx.abs() * PM.normalise();
/*  A+ku=P;
  AP*PM=0;
  ku*(M-A-ku)=0;
  k*(u*(M-A))-k2*u2=0;
  
  M-*/
}

void SurfaceProjection::updateBackgroundGridInfo_original()
{
  getAllCells(m_Cells, m_BGrid);
  getNodesFromCells(m_Cells, m_Nodes, m_BGrid);
  QVector<int> m_LNodes(m_Nodes.size());
  for (int i = 0; i < m_LNodes.size(); ++i) {
    m_LNodes[i] = i;
  }
  createNodeToNode(m_Cells, m_Nodes, m_LNodes, m_N2N, m_BGrid);
  m_EdgeLength.fill(1e99, m_BGrid->GetNumberOfPoints());
  foreach (vtkIdType id_node, m_Nodes) {
    vec3_t x;
    m_BGrid->GetPoints()->GetPoint(id_node, x.data());
    foreach (vtkIdType id_neigh, m_N2N[id_node]) {
      vec3_t xn;
      m_BGrid->GetPoints()->GetPoint(id_neigh, xn.data());
      m_EdgeLength[id_node] = min(m_EdgeLength[id_node], (x-xn).abs());
    }
    if (m_N2N[id_node].size() < 2) {
      EG_BUG;
    }
  }
  // create m_Triangles
  m_Triangles.resize(m_BGrid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    vtkIdType Npts, *pts;
    m_BGrid->GetCellPoints(id_cell, Npts, pts);
    if (Npts == 3) {
      m_BGrid->GetPoints()->GetPoint(pts[0], m_Triangles[id_cell].m_a.data());
      m_BGrid->GetPoints()->GetPoint(pts[1], m_Triangles[id_cell].m_b.data());
      m_BGrid->GetPoints()->GetPoint(pts[2], m_Triangles[id_cell].m_c.data());
      m_Triangles[id_cell].m_id_a = pts[0];
      m_Triangles[id_cell].m_id_b = pts[1];
      m_Triangles[id_cell].m_id_c = pts[2];
      m_Triangles[id_cell].m_g1 = m_Triangles[id_cell].m_b - m_Triangles[id_cell].m_a;
      m_Triangles[id_cell].m_g2 = m_Triangles[id_cell].m_c - m_Triangles[id_cell].m_a;
      m_Triangles[id_cell].m_g3 = m_Triangles[id_cell].m_g1.cross(m_Triangles[id_cell].m_g2);
      m_Triangles[id_cell].m_A  = 0.5*m_Triangles[id_cell].m_g3.abs();
      m_Triangles[id_cell].m_g3.normalise();
      m_Triangles[id_cell].m_G.column(0, m_Triangles[id_cell].m_g1);
      m_Triangles[id_cell].m_G.column(1, m_Triangles[id_cell].m_g2);
      m_Triangles[id_cell].m_G.column(2, m_Triangles[id_cell].m_g3);
      m_Triangles[id_cell].m_GI = m_Triangles[id_cell].m_G.inverse();
      m_Triangles[id_cell].m_smallest_length = (m_Triangles[id_cell].m_b - m_Triangles[id_cell].m_a).abs();
      m_Triangles[id_cell].m_smallest_length = min(m_Triangles[id_cell].m_smallest_length, (m_Triangles[id_cell].m_c - m_Triangles[id_cell].m_b).abs());
      m_Triangles[id_cell].m_smallest_length = min(m_Triangles[id_cell].m_smallest_length, (m_Triangles[id_cell].m_a - m_Triangles[id_cell].m_c).abs());
    } else {
      EG_ERR_RETURN("only triangles allowed at the moment");
    }
  }
  
  // compute node normals
  m_NodeNormals.resize(m_BGrid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node] = vec3_t(0,0,0);
  }
  foreach (Triangle T, m_Triangles) {
    m_NodeNormals[T.m_id_a] += T.m_A*T.m_g3;
    m_NodeNormals[T.m_id_b] += T.m_A*T.m_g3;
    m_NodeNormals[T.m_id_c] += T.m_A*T.m_g3;
  }
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node].normalise();
  }
  
  // compute maximum angle per node
  QVector<double> min_cos(m_BGrid->GetNumberOfPoints(), 1.0);
  foreach (Triangle T, m_Triangles) {
    double cosa = T.m_g3*m_NodeNormals[T.m_id_a];
    double cosb = T.m_g3*m_NodeNormals[T.m_id_b];
    double cosc = T.m_g3*m_NodeNormals[T.m_id_c];
    min_cos[T.m_id_a] = min(cosa, min_cos[T.m_id_a]);
    min_cos[T.m_id_b] = min(cosb, min_cos[T.m_id_b]);
    min_cos[T.m_id_c] = min(cosc, min_cos[T.m_id_c]);
  }
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    double s = sqrt(1.0 - sqr(min(1 - 1e-20, min_cos[id_node])));
    m_EdgeLength[id_node] *= m_RadiusFactor*min_cos[id_node]/s;
  }
  
}

void SurfaceProjection::updateBackgroundGridInfo()
{
  EG_VTKDCN(vtkCharArray, node_type, m_BGrid, "node_type");//node type
  
  getAllCells(m_Cells, m_BGrid);
  getNodesFromCells(m_Cells, m_Nodes, m_BGrid);
  
  setBoundaryCodes(GuiMainWindow::pointer()->getAllBoundaryCodes());
//   qDebug()<<"getBoundaryCodes()="<<getBoundaryCodes();
  
  setAllCells();
  readVMD();
  
  UpdatePotentialSnapPoints(true,false);
  l2l_t  n2n   = getPartN2N();
  l2l_t  n2c = getPartN2C();
  l2l_t  c2c   = getPartC2C();
  l2g_t cells = getPartCells();
  g2l_t _cells = getPartLocalCells();
  l2g_t nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  
//   qDebug()<<"getBoundaryCodes()="<<getBoundaryCodes();
  
  QVector<int> m_LNodes(m_Nodes.size());
  for (int i = 0; i < m_LNodes.size(); ++i) {
    m_LNodes[i] = i;
  }
  
  createNodeToNode(m_Cells, m_Nodes, m_LNodes, m_N2N, m_BGrid);
  
  m_EdgeLength.fill(1e99, m_BGrid->GetNumberOfPoints());
  foreach (vtkIdType id_node, m_Nodes) {
    vec3_t x;
    m_BGrid->GetPoints()->GetPoint(id_node, x.data());
    foreach (vtkIdType id_neigh, m_N2N[id_node]) {
      vec3_t xn;
      m_BGrid->GetPoints()->GetPoint(id_neigh, xn.data());
      m_EdgeLength[id_node] = min(m_EdgeLength[id_node], (x-xn).abs());
    }
    if (m_N2N[id_node].size() < 2) {
      EG_BUG;
    }
  }
  // create m_Triangles
  m_Triangles.resize(m_BGrid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    m_Triangles[id_cell] = Triangle(m_BGrid, id_cell);
    // TODO: Store infos about neighbour cells for each triangle
    
    for(int i=0;i<3;i++) {
      int i_cell = _cells[id_cell];
      if(c2c[i_cell][i]<0) {
        m_Triangles[id_cell].m_has_neighbour[i] = false;
      }
      else {
        m_Triangles[id_cell].m_has_neighbour[i] = true;
      }
    }
    
  }
  
  // compute node normals
  m_NodeNormals.resize(m_BGrid->GetNumberOfPoints());
//   QVector <double> NodeAngles(m_NodeNormals.size());
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node] = vec3_t(0,0,0);
//     NodeAngles[id_node]=0;
  }
  
  foreach (Triangle T, m_Triangles) {
    double angle_a = GeometryTools::angle(m_BGrid,T.m_id_c,T.m_id_a,T.m_id_b);
    double angle_b = GeometryTools::angle(m_BGrid,T.m_id_a,T.m_id_b,T.m_id_c);
    double angle_c = GeometryTools::angle(m_BGrid,T.m_id_b,T.m_id_c,T.m_id_a);
    if(isnan(angle_a) || isinf(angle_a)) EG_BUG;
    if(isnan(angle_b) || isinf(angle_b)) EG_BUG;
    if(isnan(angle_c) || isinf(angle_c)) EG_BUG;
    if(!checkVector(T.m_g3)) {
      qWarning()<<"T.m_g3="<<T.m_g3;
      EG_BUG;
    }
    double total_angle = angle_a + angle_b + angle_c;
    m_NodeNormals[T.m_id_a] += angle_a*T.m_g3;
    m_NodeNormals[T.m_id_b] += angle_b*T.m_g3;
    m_NodeNormals[T.m_id_c] += angle_c*T.m_g3;
    if(!checkVector(m_NodeNormals[T.m_id_a])) EG_BUG;
    if(!checkVector(m_NodeNormals[T.m_id_b])) EG_BUG;
    if(!checkVector(m_NodeNormals[T.m_id_c])) EG_BUG;
  }
  
//   qDebug()<<"===STARTING NORMAL CALCULATION===";
  
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
//     qDebug()<<"id_node="<<id_node<<" and node_type="<< VertexType2Str(node_type->GetValue(id_node));
//     qDebug()<<"n2n["<<id_node<<"]="<<n2n[id_node];
    
    // take into account curvature at boundaries
    if( false && node_type->GetValue(id_node)==VTK_BOUNDARY_EDGE_VERTEX) {
//       qDebug()<<"looking for edges...";
      QVector <vtkIdType> id_snappoints = getPotentialSnapPoints(id_node);
//       qDebug()<<"id_snappoints.size()="<<id_snappoints.size();
//       qDebug()<<"id_snappoints[0]="<<id_snappoints[0];
//       qDebug()<<"id_snappoints[1]="<<id_snappoints[1];
      vec3_t x0,x1,x2;
      m_Grid->GetPoints()->GetPoint(id_node, x0.data());
      m_Grid->GetPoints()->GetPoint(id_snappoints[0], x1.data());
      m_Grid->GetPoints()->GetPoint(id_snappoints[1], x2.data());
//       qDebug()<<"x0="<<x0<<" x1="<<x1<<" x="<<x2;
      
//       t1.cross(n1);
//       t2.cross(n2);
      
/*      QVector <vtkIdType> edge_cells_1;
      getEdgeCells(id_node,id_snappoints[0],edge_cells_1);
      vtkIdType id_cell = edge_cells_1[0];
      int side = getSide(id_cell,m_BGrid,id_node,id_snappoints[0]);
      vtkIdType *pts, N_pts;
      m_BGrid->GetCellPoints(id_cell, N_pts, pts);
      vec3_t t1;
      if(pts[side]==id_node) {
        t1 = x1-x0;
      }
      else {
        t1 = x0-x1;
      }
      vec3_t n1 = m_Triangles[id_cell].m_g3;
      vec3_t Nedge1 = t1.cross(n1);
      qDebug()<<"Nedge1="<<Nedge1;*/
      
//       QVector <int> i_cell_vector = n2c[_nodes[id_node]];
//       vtkIdType id_cell = cells[i_cell_vector[0]];
      
      vec3_t Nedge1 = getEdgeNormal(id_node, id_snappoints[0]);
      vec3_t Nedge2 = getEdgeNormal(id_node, id_snappoints[1]);
      vec3_t N = Nedge1+Nedge2;//(x0-x1) + (x0-x2);
      
      m_NodeNormals[id_node] = N;
//       qDebug()<<"x0="<<x0[0]<<x0[1]<<x0[2];
//       qDebug()<<"x1="<<x1[0]<<x1[1]<<x1[2];
//       qDebug()<<"x2="<<x2[0]<<x2[1]<<x2[2];
    }
    m_NodeNormals[id_node].normalise();
//     qDebug()<<"m_NodeNormals["<<id_node<<"]="<<m_NodeNormals[id_node];
  }
  
  // get the control points
  ///------------------------------
  /// UNDER CONSTRUCTION
  ///------------------------------
  
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    Triangle T = m_Triangles[id_cell];
    for(int i_edge = 0; i_edge < 3; i_edge++) {
/*      if(!m_ControlPoints.contains(OrderedPair(p1,p2))) {
      
      }*/
    }
  }
  
  for(int i_tri=0; i_tri<m_Triangles.size(); i_tri++) {
    m_Triangles[i_tri].m_Normal_a = m_NodeNormals[m_Triangles[i_tri].m_id_a];
    m_Triangles[i_tri].m_Normal_b = m_NodeNormals[m_Triangles[i_tri].m_id_b];
    m_Triangles[i_tri].m_Normal_c = m_NodeNormals[m_Triangles[i_tri].m_id_c];
    
    Triangle T = m_Triangles[i_tri];
    vec3_t X_200 = T.m_a;
    vec3_t X_020 = T.m_b;
    vec3_t X_002 = T.m_c;
    vec3_t X_011, X_101, X_110;
    getControlPoints_nonorthogonal(T,X_011, X_101, X_110, 1e99);
    m_ControlPoints[OrderedPair(T.m_id_b, T.m_id_c)] = X_011;
    m_ControlPoints[OrderedPair(T.m_id_c, T.m_id_a)] = X_101;
    m_ControlPoints[OrderedPair(T.m_id_a, T.m_id_b)] = X_110;
    
    /*
    foreach(triangle) {
      foreach(side) {
        if(side not done) {
          foreach(triangle in stencil) {
            foreach(different edge and point) {
              if(intersection) x=intersection
            }
          }
        }
      }
    }
    */
  }
  ///------------------------------
  
  // store the bezier triangles
  m_BezierTriangles.resize(m_Triangles.size());
  for(int i_tri=0; i_tri<m_Triangles.size(); i_tri++) {
    Triangle T = m_Triangles[i_tri];
    vec3_t X_200 = T.m_a;
    vec3_t X_020 = T.m_b;
    vec3_t X_002 = T.m_c;
    vec3_t X_011 = m_ControlPoints[OrderedPair(T.m_id_b, T.m_id_c)];
    vec3_t X_101 = m_ControlPoints[OrderedPair(T.m_id_c, T.m_id_a)];
    vec3_t X_110 = m_ControlPoints[OrderedPair(T.m_id_a, T.m_id_b)];
    
    m_BezierTriangles[i_tri] = BezierTriangle(X_200, X_020, X_002, X_011, X_101, X_110);
    m_BezierTriangles[i_tri].m_has_neighbour = m_Triangles[i_tri].m_has_neighbour;
  }
  
  // compute maximum angle per node
  QVector<double> min_cos(m_BGrid->GetNumberOfPoints(), 1.0);
  foreach (Triangle T, m_Triangles) {
    double cosa = T.m_g3*m_NodeNormals[T.m_id_a];
    double cosb = T.m_g3*m_NodeNormals[T.m_id_b];
    double cosc = T.m_g3*m_NodeNormals[T.m_id_c];
    min_cos[T.m_id_a] = min(cosa, min_cos[T.m_id_a]);
    min_cos[T.m_id_b] = min(cosb, min_cos[T.m_id_b]);
    min_cos[T.m_id_c] = min(cosc, min_cos[T.m_id_c]);
  }
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    double s = sqrt(1.0 - sqr(min(1 - 1e-20, min_cos[id_node])));
    m_EdgeLength[id_node] *= m_RadiusFactor*min_cos[id_node]/s;
  }
}

vec3_t SurfaceProjection::correctCurvature2(int i_tri, vec3_t g_M)
{
  return m_BezierTriangles[i_tri].projectOnQuadraticBezierTriangle(g_M);
}// end of correctCurvature2
