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

#include <vtkUnstructuredGridWriter.h>

///@@@ TODO: Delete those grids somewhere
SurfaceProjection::SurfaceProjection()
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
      double scal = (xp - T.a)*T.g3;
      vec3_t x1, x2;
      if (scal > 0) {
        x1 = xp + T.g3;
        x2 = xp - scal*T.g3 - T.g3;
      } else {
        x1 = xp - T.g3;
        x2 = xp - scal*T.g3 + T.g3;
      }
      double d = 1e99;

      bool intersects_face = GeometryTools::intersectEdgeAndTriangle(T.a, T.b, T.c, x1, x2, xi, ri);
      if (!intersects_face) {
        double kab = GeometryTools::intersection(T.a, T.b - T.a, xp, T.b - T.a);
        double kac = GeometryTools::intersection(T.a, T.c - T.a, xp, T.c - T.a);
        double kbc = GeometryTools::intersection(T.b, T.c - T.b, xp, T.c - T.b);
        double dab = (T.a + kab*(T.b-T.a) - xp).abs();
        double dac = (T.a + kac*(T.c-T.a) - xp).abs();
        double dbc = (T.b + kbc*(T.c-T.b) - xp).abs();
        bool set = false;
        if ((kab >= 0) && (kab <= 1)) {
          if (dab < d) {
            xi = T.a + kab*(T.b-T.a);
            d = dab;
            set = true;
          }
        }
        if ((kac >= 0) && (kac <= 1)) {
          if (dac < d) {
            xi = T.a + kac*(T.c-T.a);
            d = dac;
            set = true;
          }
        }
        if ((kbc >= 0) && (kbc <= 1)) {
          if (dbc < d) {
            xi = T.b + kbc*(T.c-T.b);
            d = dbc;
            set = true;
          }
        }
        double da = (T.a - xp).abs();
        double db = (T.b - xp).abs();
        double dc = (T.c - xp).abs();
        if (da < d) {
          xi = T.a;
          d = da;
          set = true;
        }
        if (db < d) {
          xi = T.b;
          d = db;
        }
        if (dc < d) {
          xi = T.c;
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
      double g = dx*T.g3;
      if (intersects_face) {
        d = fabs(g);
      }
      double w = 1;
      w *= m_DistWeight*pow(L/max(1e-6*L, d), m_DistExp);
      if (dx.abs() < 1e-6*L) {
        w *= m_DirWeight;
      } else {
        dx.normalise();
        w *= m_DirWeight*pow(fabs(dx*T.g3), m_DirExp);
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

// #define EGVTKOBJECT_CREATENODEFIELD(FIELD,TYPE,OW)
// if (!grid->GetPointData()->GetArray(FIELD)) {
// EG_VTKSP(TYPE, var);
// var->SetName(FIELD);
// var->SetNumberOfValues(Nnodes);
// grid->GetPointData()->AddArray(var);
// for (int i = 0; i < grid->GetNumberOfPoints(); ++i) {
// var->SetValue(i,0);
// }
// } else if (OW) {
// EG_VTKDCN(TYPE, var, grid, FIELD);
// var->SetNumberOfValues(Nnodes);
// for (int i = 0; i < grid->GetNumberOfPoints(); ++i) {
// var->SetValue(i,0);
// }
// }

void SurfaceProjection::writeGridWithNormals()
{
  //qDebug()<<"void SurfaceProjection::writeGridWithNormals() called";
  
  vtkDoubleArray *vectors = vtkDoubleArray::New();
  vectors->SetName("normals");
  vectors->SetNumberOfComponents(3);
  vectors->SetNumberOfTuples(m_BGrid->GetNumberOfPoints());
  
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    vec3_t N = m_NodeNormals[id_node];
    double n[3];
    n[0]=N[0];
    n[1]=N[1];
    n[2]=N[2];
    vectors->InsertTuple(id_node,n);
  }
  
  m_BGrid->GetPointData()->SetVectors(vectors);
/*  vectors->SetName("normals2");
  m_BGrid->GetPointData()->SetVectors(vectors);
  vectors->SetName("normals3");
  m_BGrid->GetPointData()->SetVectors(vectors);*/
//   m_BGrid->GetPointData()->AddArray(vectors);
  
  vectors->Delete();
  
  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
  vtu->SetFileName("m_BGrid.vtu");
  vtu->SetDataModeToBinary();
  vtu->SetInput(m_BGrid);
  vtu->Write();
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

double interpolate(vec2_t A, vec2_t nA, vec2_t M, vec2_t I, vec2_t nI)
{
  bool DEBUG = false;
  
  if(DEBUG) qDebug()<<"double interpolate(vec2_t A, vec2_t nA, vec2_t M, vec2_t I, vec2_t nI) called";
  if(DEBUG) cout<<"A="<<A<<endl;
  if(DEBUG) cout<<"nA="<<nA<<endl;
  if(DEBUG) cout<<"M="<<M<<endl;
  if(DEBUG) cout<<"I="<<I<<endl;
  if(DEBUG) cout<<"nI="<<nI<<endl;
  
  double ret = 0;
  
  double alpha0 = -nA[0]/nA[1];
  double alpha1 = -nI[0]/nI[1];
  
  // f(x)=a*x^3 + b*x^2 + c*x + d
  double a,b,c,d;
  d = 0;
  c = alpha0;
  b = -((alpha1-c)-3*(-c));
  a = -b-c;
  
  double xM = M[0];
  
  if(DEBUG) cout<<"a="<<a<<endl;
  if(DEBUG) cout<<"b="<<b<<endl;
  if(DEBUG) cout<<"c="<<c<<endl;
  if(DEBUG) cout<<"d="<<d<<endl;
  ret = a*pow(xM,3) + b*pow(xM,2) + c*xM + d;
  
  // B(t)=(1-t^3)*P0 + 3*(1-t)^2*t*P1 + 3*(1-t)*t^2*P2 + t^3*P3;
  
//   return -(xM*xM) + xM;
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

vec2_t projectVectorOnPlane(vec3_t V,vec3_t i,vec3_t j)
{
  double x = V*i/i.abs2();
  double y = V*j/j.abs2();
  return vec2_t(x,y);
}

vec3_t SurfaceProjection::correctCurvature(int i_tri, vec3_t r)
{
  // initialization
  bool DEBUG = false;
  if(DEBUG) qDebug()<<"vec3_t SurfaceProjection::correctCurvature(int i_tri, vec3_t r) called";
  vec3_t x(0,0,0);
  
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
  
  // before knowing intersections
  
  Triangle T = m_Triangles[i_tri];
  vec3_t g_A = T.a;
  vec3_t g_B = T.b;
  vec3_t g_C = T.c;
  vec3_t g_M = g_A+T.G*r;
  
  vec3_t l_A(0,0,0);
  vec3_t l_B(1,0,0);
  vec3_t l_C(0,1,0);
  vec3_t l_M = r;
  
  vec2_t t_A(0,0);
  vec2_t t_B(1,0);
  vec2_t t_C(0,1);
  vec2_t t_M(r[0],r[1]);
  
  vec3_t g_nA = m_NodeNormals[T.id_a];
  vec3_t g_nB = m_NodeNormals[T.id_b];
  vec3_t g_nC = m_NodeNormals[T.id_c];
  
  vec2_t pm1_A(0,0);
  vec2_t pm1_I1(1,0);
//   vec2_t pm1_(0,0);
  
  vec2_t pm2_B(0,0);
  vec2_t pm2_I2(1,0);
  
  vec2_t pm3_C(0,0);
  vec2_t pm3_I3(1,0);
  
  vec3_t g_g1 = T.g1;
  vec3_t g_g2 = T.g2;
  vec3_t g_g3 = T.g3;
  
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
  
  vec3_t g_I1 = g_A+T.G*l_I1;
  vec3_t g_I2 = g_A+T.G*l_I2;
  vec3_t g_I3 = g_A+T.G*l_I3;
  
  vec3_t tmp;
  tmp = (g_nI1); vec3_t l_nI1 = T.GI*tmp;
  tmp = (g_nI2); vec3_t l_nI2 = T.GI*tmp;
  tmp = (g_nI3); vec3_t l_nI3 = T.GI*tmp;
  
  tmp = (g_nA); vec3_t l_nA = T.GI*tmp;
  tmp = (g_nB); vec3_t l_nB = T.GI*tmp;
  tmp = (g_nC); vec3_t l_nC = T.GI*tmp;
  
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
  
  // plane 1: BC, nI3 -> pnoe1_
  // plane 2: CA, nI1 -> pnoe2_
  // plane 3: AB, nI2 -> pnoe3_
  
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
  
  vec2_t pnoe1_K1;
  if(!intersection(k1,k2,pnoe1_B,pnoe1_tB,pnoe1_C,pnoe1_tC)) {
    //cout<<"pnoe1_B="<<pnoe1_B<<endl;
    //cout<<"pnoe1_tB="<<pnoe1_tB<<endl;
    //cout<<"pnoe1_C="<<pnoe1_C<<endl;
    //cout<<"pnoe1_tC="<<pnoe1_tC<<endl;
//     EG_BUG;
    pnoe1_K1 = 0.5*(pnoe1_B + pnoe1_C);
  }
  else {
    pnoe1_K1 = pnoe1_B + k1*pnoe1_tB;
  }
  
  vec3_t l_K1 = l_B + pnoe1_K1[0]*l_BC + pnoe1_K1[1]*l_nI1;
  ///////////////
  vec2_t pnoe2_C(0,0);
  vec2_t pnoe2_A(1,0);
  vec2_t pnoe2_nC = projectVectorOnPlane(l_nC,l_CA,l_nI2);
  vec2_t pnoe2_nA = projectVectorOnPlane(l_nA,l_CA,l_nI2);
  vec2_t pnoe2_nI2(0,1);
  
  vec2_t pnoe2_tC = turnRight(pnoe2_nC);
  vec2_t pnoe2_tA = turnRight(pnoe2_nA);
  
  vec2_t pnoe2_K2;
  if(!intersection(k1,k2,pnoe2_C,pnoe2_tC,pnoe2_A,pnoe2_tA)) {
    //cout<<"pnoe2_C="<<pnoe2_C<<endl;
    //cout<<"pnoe2_tC="<<pnoe2_tC<<endl;
    //cout<<"pnoe2_A="<<pnoe2_A<<endl;
    //cout<<"pnoe2_tA="<<pnoe2_tA<<endl;
//     EG_BUG;
    pnoe2_K2 = 0.5*(pnoe2_C + pnoe2_A);
  }
  else {
    pnoe2_K2 = pnoe2_C + k1*pnoe2_tC;
  }
  
  vec3_t l_K2 = l_C + pnoe2_K2[0]*l_CA + pnoe2_K2[1]*l_nI2;
  ///////////////
  vec2_t pnoe3_A(0,0);
  vec2_t pnoe3_B(1,0);
  vec2_t pnoe3_nA = projectVectorOnPlane(l_nA,l_AB,l_nI3);
  vec2_t pnoe3_nB = projectVectorOnPlane(l_nB,l_AB,l_nI3);
  vec2_t pnoe3_nI3(0,1);
  
  vec2_t pnoe3_tA = turnRight(pnoe3_nA);
  vec2_t pnoe3_tB = turnRight(pnoe3_nB);
  
  vec2_t pnoe3_K3;
  if(!intersection(k1,k2,pnoe3_A,pnoe3_tA,pnoe3_B,pnoe3_tB)) {
    //cout<<"pnoe3_A="<<pnoe3_A<<endl;
    //cout<<"pnoe3_tA="<<pnoe3_tA<<endl;
    //cout<<"pnoe3_B="<<pnoe3_B<<endl;
    //cout<<"pnoe3_tB="<<pnoe3_tB<<endl;
//     EG_BUG;
    pnoe3_K3 = 0.5*(pnoe3_A + pnoe3_B);
  }
  else {
    pnoe3_K3 = pnoe3_A + k1*pnoe3_tA;
  }
  
  
  vec3_t l_K3 = l_A + pnoe3_K3[0]*l_AB + pnoe3_K3[1]*l_nI3;
  ///////////////
  
  vec3_t g_K1 = g_A+T.G*l_K1;
  vec3_t g_K2 = g_A+T.G*l_K2;
  vec3_t g_K3 = g_A+T.G*l_K3;
  
  vec3_t g_J1;
  vec3_t g_J2;
  vec3_t g_J3;
  
  //qDebug()<<"=== ORTHOGONAL PLANES ===";
  getControlPoints_orthogonal(T,g_J1,g_J2,g_J3);
  //qDebug()<<"=== NON-ORTHOGONAL PLANES ===";
  getControlPoints_nonorthogonal(T,g_K1,g_K2,g_K3);

  vec3_t X_200 = g_A;
  vec3_t X_020 = g_B;
  vec3_t X_002 = g_C;
  vec3_t X_011 = g_K1;
  vec3_t X_101 = g_K2;
  vec3_t X_110 = g_K3;
//   writeBezierSurface(X_200, X_020, X_002, X_011, X_101, X_110);
  

  BezierTriangle bezier_triangle(X_200, X_020, X_002, X_011, X_101, X_110);
  return bezier_triangle.projectOnQuadraticBezierTriangle(g_M);
//   return bezier_triangle.QuadraticBezierTriangle(t_M);
  
/*  intersection(k1,k2,pm1_A,pm1_nA,pm1_B,pm1_nB);
  intersection(k1,k2,pm2_B,pm2_nB,pm2_C,pm2_nC);
  intersection(k1,k2,pm3_C,pm3_nC,pm3_A,pm3_nA);*/
  
  // interpolation attempts
  
/*  double z1 = interpolate(pm1_A, pm1_nA, pm1_M, pm1_I1, pm1_nI1);
  double z2 = interpolate(pm2_B, pm2_nB, pm2_M, pm2_I2, pm2_nI2);
  double z3 = interpolate(pm3_C, pm3_nC, pm3_M, pm3_I3, pm3_nI3);
  double z = (z1+z2+z3)/3.0;*/
  
  vec2_t pm1_Z1 = interpolate_bezier_normals(pm1_M[0],pm1_A,pm1_nA,pm1_I1,pm1_nI1);
  vec2_t pm2_Z2 = interpolate_bezier_normals(pm2_M[0],pm2_B,pm2_nB,pm2_I2,pm2_nI2);
  vec2_t pm3_Z3 = interpolate_bezier_normals(pm3_M[0],pm3_C,pm3_nC,pm3_I3,pm3_nI3);
  
  double z1 = pm1_Z1[1];
  double z2 = pm2_Z2[1];
  double z3 = pm3_Z3[1];
  double z = z1;//(z1+z2+z3)/3.0;
  
/*  vec3_t g_Z1 = interpolate_bezier_normals(pm1_M[0],g_A,g_nA,g_I1,g_nI1,g_AI1,g_g3);
  vec3_t g_Z2 = interpolate_bezier_normals(pm2_M[0],g_B,g_nB,g_I2,g_nI2,g_BI2,g_g3);
  vec3_t g_Z3 = interpolate_bezier_normals(pm3_M[0],g_C,g_nC,g_I3,g_nI3,g_CI3,g_g3);
  return (1.0/3.0)*(g_Z1+g_Z2+g_Z3);*/
  
  vec3_t l_X = l_M + z*l_g3;
  vec3_t g_X = g_A+T.G*l_X;
  
  
  vec3_t A,M,I;
  vec3_t nA,nM,nI;
  
  A = vec3_t(pm2_B[0],pm2_B[1],0);
  M = vec3_t(pm2_M[0],pm2_M[1],0);
  I = vec3_t(pm2_I2[0],pm2_I2[1],0);
  nA = vec3_t(pm2_nB[0],pm2_nB[1],0);
  nM = vec3_t(0,1,0);
  nI = vec3_t(pm2_nI2[0],pm2_nI2[1],0);
  
  // debugging stuff
  
  if(DEBUG) cout<<"r="<<r<<endl;
  if(DEBUG) cout<<"l_g1"<<l_g1<<endl;
  if(DEBUG) cout<<"l_g2"<<l_g2<<endl;
  if(DEBUG) cout<<"l_g3"<<l_g3<<endl;
  if(DEBUG) cout<<"g_g1"<<g_g1<<endl;
  if(DEBUG) cout<<"g_g2"<<g_g2<<endl;
  if(DEBUG) cout<<"g_g3"<<g_g3<<endl;
  if(DEBUG) cout<<"l_nI1"<<l_nI1<<endl;
  if(DEBUG) cout<<"l_nI2"<<l_nI2<<endl;
  if(DEBUG) cout<<"l_nI3"<<l_nI3<<endl;
  if(DEBUG) cout<<"l_nA"<<l_nA<<endl;
  if(DEBUG) cout<<"l_nB"<<l_nB<<endl;
  if(DEBUG) cout<<"l_nC"<<l_nC<<endl;
  if(DEBUG) cout<<"l_nI2="<<l_nI2<<endl;
  if(DEBUG) cout<<"l_BI2="<<l_BI2<<endl;
  if(DEBUG) cout<<"g_nA"<<g_nA<<endl;
  if(DEBUG) cout<<"g_nI1"<<g_nI1<<endl;
  if(DEBUG) cout<<"g_nB"<<g_nB<<endl;
  if(DEBUG) cout<<"g_nI2"<<g_nI2<<endl;
  if(DEBUG) cout<<"g_nC"<<g_nC<<endl;
  if(DEBUG) cout<<"g_nI3"<<g_nI3<<endl;
  
  QVector < QPair<vec3_t,vec3_t> > points;
  points.push_back(QPair<vec3_t,vec3_t>(g_A,g_nA));
  points.push_back(QPair<vec3_t,vec3_t>(g_B,g_nB));
  points.push_back(QPair<vec3_t,vec3_t>(g_C,g_nC));
  points.push_back(QPair<vec3_t,vec3_t>(g_I1,g_nI1));
  points.push_back(QPair<vec3_t,vec3_t>(g_I2,g_nI2));
  points.push_back(QPair<vec3_t,vec3_t>(g_I3,g_nI3));
  points.push_back(QPair<vec3_t,vec3_t>(g_M,g_g3));
  
/*  points.push_back(QPair<vec3_t,vec3_t>(A,nA));
  points.push_back(QPair<vec3_t,vec3_t>(M,nM));
  points.push_back(QPair<vec3_t,vec3_t>(I,nI));*/
  
  QVector < QPair<vec3_t,vec3_t> > lines;
  lines.push_back(QPair<vec3_t,vec3_t>(g_A,g_I1));
  lines.push_back(QPair<vec3_t,vec3_t>(g_B,g_I2));
  lines.push_back(QPair<vec3_t,vec3_t>(g_C,g_I3));
  
  if(DEBUG) debug_output(points, lines);
  
  // returning value
  x = g_X;
  return x;
}

bool SurfaceProjection::projectOnTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d, const Triangle& T)
{
  xi = vec3_t(1e99,1e99,1e99);
  double scal = (xp - T.a)*T.g3;
  vec3_t x1, x2;
  if (scal > 0) {
    x1 = xp + T.g3;
    x2 = xp - scal*T.g3 - T.g3;
  } else {
    x1 = xp - T.g3;
    x2 = xp - scal*T.g3 + T.g3;
  }
  d = 1e99;
  bool intersects_face = GeometryTools::intersectEdgeAndTriangle(T.a, T.b, T.c, x1, x2, xi, ri);
  if (intersects_face) {
    vec3_t dx = xp - T.a;
    d = fabs(dx*T.g3);
  } else {
    double kab = GeometryTools::intersection(T.a, T.b - T.a, xp, T.b - T.a);
    double kac = GeometryTools::intersection(T.a, T.c - T.a, xp, T.c - T.a);
    double kbc = GeometryTools::intersection(T.b, T.c - T.b, xp, T.c - T.b);
    double dab = (T.a + kab*(T.b-T.a) - xp).abs();
    double dac = (T.a + kac*(T.c-T.a) - xp).abs();
    double dbc = (T.b + kbc*(T.c-T.b) - xp).abs();
    bool set = false;
    if ((kab >= 0) && (kab <= 1)) {
      if (dab < d) {
        xi = T.a + kab*(T.b-T.a);
        d = dab;
        set = true;
      }
    }
    if ((kac >= 0) && (kac <= 1)) {
      if (dac < d) {
        xi = T.a + kac*(T.c-T.a);
        d = dac;
        set = true;
      }
    }
    if ((kbc >= 0) && (kbc <= 1)) {
      if (dbc < d) {
        xi = T.b + kbc*(T.c-T.b);
        d = dbc;
        set = true;
      }
    }
    double da = (T.a - xp).abs();
    double db = (T.b - xp).abs();
    double dc = (T.c - xp).abs();
    if (da < d) {
      xi = T.a;
      d = da;
      set = true;
    }
    if (db < d) {
      xi = T.b;
      d = db;
    }
    if (dc < d) {
      xi = T.c;
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
  return intersects_face;
}

vec3_t SurfaceProjection::cylinder(vec3_t center, double radius, int i_tri, vec3_t r)
{
  Triangle T = m_Triangles[i_tri];
  vec3_t g_A = T.a;
  vec3_t g_M = g_A+T.G*r;
  vec3_t P;
  P[2]=g_M[2];
  double L = g_M[0]*g_M[0]+g_M[1]*g_M[1];
  P[0]=radius * g_M[0]/L;
  P[1]=radius * g_M[1]/L;
  return P;
}

vec3_t SurfaceProjection::projectWithGeometry(vec3_t xp, vtkIdType id_node)
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
      bool intersects = projectOnTriangle(xp, xi, ri, d, T);
      if (!intersects || (d > 0.1*T.smallest_length)) {
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
      bool intersects = projectOnTriangle(xp, xi, ri, d, T);
      if (d < d_min) {
        x_proj = xi;
        r_proj = ri;
        d_min = d;
        m_ProjTriangles[id_node] = i_triangles;
        on_triangle = intersects;
      }
    }
    if (x_proj[0] > 1e98) {
      writeGrid(GuiMainWindow::pointer()->getGrid(),"griddump");
      EG_BUG;
    }
  }
  if(on_triangle) {
    vec3_t center(0,0,0);
    double radius = 1;
    x_proj = cylinder(center, radius, m_ProjTriangles[id_node], r_proj);
//     x_proj = correctCurvature(m_ProjTriangles[id_node], r_proj);
  }
  return x_proj;
}

vec3_t SurfaceProjection::project(vec3_t x, vtkIdType id_node)
{
  if (m_UseLevelSet) {
    x = projectWithLevelSet(x);
  } else {
    if (id_node < 0) {
      EG_BUG;
    }
    x = projectWithGeometry(x, id_node);
  }
//   writeGrid(m_FGrid,"m_FGrid");
//   writeGrid(m_BGrid,"m_BGrid");
  return x;
}

vec3_t intersectionOnPlane(vec3_t v, vec3_t A, vec3_t nA, vec3_t B, vec3_t nB)
{
  vec3_t u = B-A;
//   u.normalise();
  v.normalise();
  v = u.abs()*v;
  
  //cout<<"u="<<u<<" v="<<v<<endl;

  vec2_t p_A(0,0);
  vec2_t p_B(1,0);
  vec2_t p_nA = projectVectorOnPlane(nA,u,v);
  vec2_t p_nB = projectVectorOnPlane(nB,u,v);
  
  vec2_t p_tA = turnRight(p_nA);
  vec2_t p_tB = turnRight(p_nB);
  
  double k1, k2;
  vec2_t p_K;
  if(!intersection(k1, k2, p_A, p_tA, p_B, p_tB)) {
    //qDebug()<<"WARNING: No intersection found!!!";
    p_K = 0.5*(p_A + p_B);
  }
  else {
    p_K = p_A + k1*p_tA;
  }
  
  //cout<<"nA="<<nA<<endl;
  //cout<<"p_nA="<<p_nA<<endl;
  //cout<<"p_tA="<<p_tA<<endl;
  //cout<<"p_K="<<p_K<<endl;
  vec3_t K = A + p_K[0]*u + p_K[1]*v;
  //cout<<"K="<<K<<endl;
  return K;
}

int SurfaceProjection::getControlPoints_orthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110)
{
  vec3_t A=T.a;
  vec3_t B=T.b;
  vec3_t C=T.c;
  vec3_t nA = m_NodeNormals[T.id_a];
  vec3_t nB = m_NodeNormals[T.id_b];
  vec3_t nC = m_NodeNormals[T.id_c];
  
  //cout<<"nA="<<nA<<endl;
  //cout<<"nB="<<nB<<endl;
  //cout<<"nC="<<nC<<endl;

//   cout<<"A="<<A<<" B="<<B<<" C="<<C<<endl;
  //cout<<"-->BC"<<endl;
  X_011 = intersectionOnPlane(T.g3, B, nB, C, nC);
  //cout<<"-->CA"<<endl;
  X_101 = intersectionOnPlane(T.g3, C, nC, A, nA);
  //cout<<"-->AB"<<endl;
  X_110 = intersectionOnPlane(T.g3, A, nA, B, nB);
  return(0);
}

int SurfaceProjection::getControlPoints_nonorthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110)
{
  vec3_t A=T.a;
  vec3_t B=T.b;
  vec3_t C=T.c;
  vec3_t nA = m_NodeNormals[T.id_a];
  vec3_t nB = m_NodeNormals[T.id_b];
  vec3_t nC = m_NodeNormals[T.id_c];
  
//   cout<<"A="<<A<<" B="<<B<<" C="<<C<<endl;
  X_011 = intersectionOnPlane(0.5*(nB+nC), B, nB, C, nC);
  X_101 = intersectionOnPlane(0.5*(nC+nA), C, nC, A, nA);
  X_110 = intersectionOnPlane(0.5*(nA+nB), A, nA, B, nB);
  return(0);
}

void SurfaceProjection::setupInterpolationGrid()
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
  
  allocateGrid(m_BezierGrid, m_Triangles.size()*N_cells_per_triangle, m_Triangles.size()*N_points_per_triangle);

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
    getControlPoints_orthogonal(T,J1,J2,J3);
    //qDebug()<<"=== NON-ORTHOGONAL PLANES ===";
    getControlPoints_nonorthogonal(T,K1,K2,K3);
    
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
      //cout<<"A="<<T.a<<" B="<<T.b<<" C="<<T.c<<endl;
      //cout<<"J1="<<J1<<" K1="<<K1<<endl;
      //cout<<"J2="<<J2<<" K2="<<K2<<endl;
      //cout<<"J3="<<J3<<" K3="<<K3<<endl;
      //cout<<"+++++++++++++++++++++++++"<<endl;
    }
    
    BezierTriangle bezier_triangle(T.a, T.b, T.c, K1, K2, K3);
    offset += bezier_triangle.addBezierSurface(m_BezierGrid, offset, N);

    vtkIdType polyline_ortho[7];
    vtkIdType polyline_nonortho[7];
    
    polyline_ortho[0]=T.id_a;
    polyline_ortho[1]=idx_J3;
    polyline_ortho[2]=T.id_b;
    polyline_ortho[3]=idx_J1;
    polyline_ortho[4]=T.id_c;
    polyline_ortho[5]=idx_J2;
    polyline_ortho[6]=T.id_a;
    
    polyline_nonortho[0]=T.id_a;
    polyline_nonortho[1]=idx_K3;
    polyline_nonortho[2]=T.id_b;
    polyline_nonortho[3]=idx_K1;
    polyline_nonortho[4]=T.id_c;
    polyline_nonortho[5]=idx_K2;
    polyline_nonortho[6]=T.id_a;
    
    m_InterpolationGrid->InsertNextCell(4,7,polyline_ortho);cell_count++;
    m_InterpolationGrid->InsertNextCell(4,7,polyline_nonortho);cell_count++;
    
  }
  
  //qDebug()<<"node_count="<<node_count;
  //qDebug()<<"cell_count="<<cell_count;
  //qDebug()<<"offset="<<offset;
  
  EG_VTKSP(vtkUnstructuredGridWriter,vtu1);
  vtu1->SetFileName("m_InterpolationGrid.vtk");
  vtu1->SetInput(m_InterpolationGrid);
  vtu1->Write();

  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu2);
  vtu2->SetFileName("m_InterpolationGrid.vtu");
  vtu2->SetDataModeToBinary();
//   vtu2->SetDataModeToAscii();
  vtu2->SetInput(m_InterpolationGrid);
  vtu2->Write();

  EG_VTKSP(vtkUnstructuredGridWriter,vtu3);
  vtu3->SetFileName("m_BezierGrid.vtk");
  vtu3->SetInput(m_BezierGrid);
  vtu3->Write();

  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu4);
  vtu4->SetFileName("m_BezierGrid.vtu");
  vtu4->SetDataModeToBinary();
//   vtu4->SetDataModeToAscii();
  vtu4->SetInput(m_BezierGrid);
  vtu4->Write();
  
  writeGrid(m_BGrid,"m_BGrid");
  
}

// mat2_t SurfaceProjection::Jacobian_Matrix(double x, double y)
// {
//   mat2_t J;
//   return J;
// }

// vec2_t BezierProjectionFunction(double x, double y)
// {
//   vec2_t F;
//   QuadraticBezierTriangle(double u, double v, double w, vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110);
//   return F;
// }
