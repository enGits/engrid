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
#include "surfaceoperation.h"

#include <vtkCharArray.h>
#include <vtkMath.h>
#include <vtkCellArray.h>
#include <vtkPolygon.h>

#include "geometrytools.h"
using namespace GeometryTools;

SurfaceOperation::SurfaceOperation()
 : Operation()
{
  m_CellLocator = NULL;
  m_ProjectionSurface = NULL;
  
  //default values for determining node types and for smoothing operations
  Convergence=0;
  NumberOfIterations=20;
  RelaxationFactor=0.01;
  FeatureEdgeSmoothing=1;//0 by default in VTK, but we need 1 to avoid the "potatoe effect" ^^
  FeatureAngle=45;
  EdgeAngle=15;
  BoundarySmoothing=1;
}

void SurfaceOperation::operate()
{

}

ostream& operator<<(ostream &out, stencil_t S)
{
  out<<"S.id_cell1="<<S.id_cell1<<" ";
  out<<"S.id_cell2="<<S.id_cell2<<" ";
  out<<"S.sameBC="<<S.sameBC<<" ";
  out<<"S.twocells="<<S.twocells<<" ";
  out<<"S.neighbour_type="<<S.neighbour_type<<" ";
  out<<"[";
  for(int i=0;i<4;i++){
    out<<S.p[i];
    if(i!=3) out<<",";
  }
  out<<"]";
  return(out);
}

stencil_t SurfaceOperation::getStencil(vtkIdType id_cell1, int j1)
{
  if(grid->GetCellType(id_cell1)!=VTK_TRIANGLE) {
    cout<<"CELL IS NOT A TRIANGLE"<<endl;
    EG_BUG;
  }
  
  //return variable
  stencil_t S;
  
  //default values:
  S.sameBC = false;
  S.twocells = false;
  S.neighbour_type = -1;
  
  //initialize first cell
  S.id_cell1 = id_cell1;
  vtkIdType N1, *pts1;
  grid->GetCellPoints(S.id_cell1, N1, pts1);
  //place points 0,1,3
  if      (j1 == 0) { S.p[0] = pts1[2]; S.p[1] = pts1[0]; S.p[3] = pts1[1]; }
  else if (j1 == 1) { S.p[0] = pts1[0]; S.p[1] = pts1[1]; S.p[3] = pts1[2]; }
  else if (j1 == 2) { S.p[0] = pts1[1]; S.p[1] = pts1[2]; S.p[3] = pts1[0]; };
  
  //initialize second cell
  S.id_cell2 = -1;
  S.p[2]=-1;
  
  //twocells
  if (c2c[_cells[id_cell1]][j1] != -1) {//if neighbour cell
    
    //twocells
    S.twocells = true;
    S.id_cell2 = cells[c2c[_cells[id_cell1]][j1]];
    
    //sameBC
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    if(cell_code->GetValue(S.id_cell1)==cell_code->GetValue(S.id_cell2)) S.sameBC = true;
    
    //neighbour_type
    S.neighbour_type = grid->GetCellType(S.id_cell2);
    if ( S.neighbour_type == VTK_TRIANGLE) {//if neighbour cell is a triangle
      vtkIdType N2, *pts2;
      grid->GetCellPoints(S.id_cell2, N2, pts2);
      
      //place point 2
      bool p2 = false;
      if (c2c[_cells[S.id_cell2]][0] != -1) {
        if (cells[c2c[_cells[S.id_cell2]][0]] == S.id_cell1) {
          S.p[2] = pts2[2];
          p2 = true;
        }
      }
      if (c2c[_cells[S.id_cell2]][1] != -1) {
        if (cells[c2c[_cells[S.id_cell2]][1]] == S.id_cell1) {
          S.p[2] = pts2[0];
          p2 = true;
        }
      }
      if (c2c[_cells[S.id_cell2]][2] != -1) {
        if (cells[c2c[_cells[S.id_cell2]][2]] == S.id_cell1) {
          S.p[2] = pts2[1];
          p2 = true;
        }
      }
      
      if (!p2) {//failed to place point 2, appears when cell1 is linked to cell2, but cell2 not to cell1
        cout<<"S.id_cell1="<<S.id_cell1<<endl;
        cout<<"S.id_cell2="<<S.id_cell2<<endl;
        createNodeToCell(cells, nodes, _nodes, n2c, grid);
        EG_BUG;
      }
    }
  }//end of if neighbour cell
  return S;
}

int SurfaceOperation::UpdateCurrentMeshDensity()
{
  if(DebugLevel>0) cout<<"===UpdateMeshDensity START==="<<endl;
  
  getAllSurfaceCells(cells,grid);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  getNodesFromCells(cells, nodes, grid);
  setGrid(grid);
  setCells(cells);
  
  if(DebugLevel>5) cout<<"cells.size()="<<cells.size()<<endl;
  
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, grid, "node_meshdensity_current");
  foreach(vtkIdType node,nodes)
  {
    node_meshdensity_current->SetValue(node, CurrentMeshDensity(node));
  }
  if(DebugLevel>0) cout<<"===UpdateMeshDensity END==="<<endl;
  return(0);
}

int SurfaceOperation::UpdateNodeType()
{
  cout<<"=== UpdateNodeType START ==="<<endl;
  //prepare
  setAllSurfaceCells();
  
  m_PotentialSnapPoints.resize(nodes.size());
  
  //initialize default values
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  foreach(vtkIdType id_node, nodes) {
    node_type->SetValue(id_node, VTK_SIMPLE_VERTEX);
    m_PotentialSnapPoints[id_node].clear();
  }
  
  cout<<"===pre-processing==="<<endl;
  int N_edges=0;
  //We loop through edges
  foreach(vtkIdType id_cell, cells) {
    vtkIdType *pts, Npts;
    grid->GetCellPoints(id_cell, Npts, pts);
    for(int i=0;i<Npts;i++) {
      
      int i_neighbour_cell = c2c[_cells[id_cell]][i];
      if( i_neighbour_cell>=0 && cells[i_neighbour_cell] < id_cell ) continue;//already visited edge
      N_edges++;
      
      vtkIdType id_node1 = pts[i];
      vtkIdType id_node2 = pts[(i+1)%Npts];
      
      //-----------------------
      //determine edge type
      char edge = getEdgeType(id_node2,id_node1);
      //-----------------------
      //determine node type pre-processing (count nb of complex edges if the node is complex, otherwise, just count the nb of edges)
      if ( edge && node_type->GetValue(id_node1) == VTK_SIMPLE_VERTEX )
      {
        m_PotentialSnapPoints[id_node1].clear();
        m_PotentialSnapPoints[id_node1].push_back(id_node2);
        node_type->SetValue(id_node1,edge);
      }
      else if ( (edge && node_type->GetValue(id_node1) == VTK_BOUNDARY_EDGE_VERTEX) ||
                (edge && node_type->GetValue(id_node1) == VTK_FEATURE_EDGE_VERTEX) ||
                (!edge && node_type->GetValue(id_node1) == VTK_SIMPLE_VERTEX ) )
      {
        m_PotentialSnapPoints[id_node1].push_back(id_node2);
        if ( node_type->GetValue(id_node1) && edge == VTK_BOUNDARY_EDGE_VERTEX )
        {
          node_type->SetValue(id_node1, VTK_BOUNDARY_EDGE_VERTEX);//VTK_BOUNDARY_EDGE_VERTEX has priority over VTK_FEATURE_EDGE_VERTEX
        }
      }
      
      if ( edge && node_type->GetValue(id_node2) == VTK_SIMPLE_VERTEX )
      {
        m_PotentialSnapPoints[id_node2].clear();
        m_PotentialSnapPoints[id_node2].push_back(id_node1);
        node_type->SetValue(id_node2, edge);
      }
      else if ( (edge && node_type->GetValue(id_node2) == VTK_BOUNDARY_EDGE_VERTEX ) ||
                (edge && node_type->GetValue(id_node2) == VTK_FEATURE_EDGE_VERTEX) ||
                (!edge && node_type->GetValue(id_node2) == VTK_SIMPLE_VERTEX ) )
      {
        m_PotentialSnapPoints[id_node2].push_back(id_node1);
        if ( node_type->GetValue(id_node2) && edge == VTK_BOUNDARY_EDGE_VERTEX )
        {
          node_type->SetValue(id_node2, VTK_BOUNDARY_EDGE_VERTEX);//VTK_BOUNDARY_EDGE_VERTEX has priority over VTK_FEATURE_EDGE_VERTEX
        }
      }
    }
  }
  
  cout<<"N_edges="<<N_edges<<endl;
  
  //-----------------------
  //determine node type post-processing
  double CosEdgeAngle = cos((double) vtkMath::RadiansFromDegrees(this->EdgeAngle));
  cout<<"===post-processing==="<<endl;
  //This time, we loop through nodes
  foreach(vtkIdType id_node, nodes) {
    if ( node_type->GetValue(id_node) == VTK_FEATURE_EDGE_VERTEX || node_type->GetValue(id_node) == VTK_BOUNDARY_EDGE_VERTEX )
    { //see how many edges; if two, what the angle is
      
      if ( !this->BoundarySmoothing && 
           node_type->GetValue(id_node) == VTK_BOUNDARY_EDGE_VERTEX )
      {
        node_type->SetValue(id_node, VTK_FIXED_VERTEX);
      }
      else if ( m_PotentialSnapPoints[id_node].size() != 2 )
      {
        node_type->SetValue(id_node, VTK_FIXED_VERTEX);
      }
      else //check angle between edges
      {
        double x1[3], x2[3], x3[3], l1[3], l2[3];
        grid->GetPoint(m_PotentialSnapPoints[id_node][0],x1);
        grid->GetPoint(id_node,x2);
        grid->GetPoint(m_PotentialSnapPoints[id_node][1],x3);
        for (int k=0; k<3; k++)
        {
          l1[k] = x2[k] - x1[k];
          l2[k] = x3[k] - x2[k];
        }
        if ( vtkMath::Normalize(l1) >= 0.0 &&
             vtkMath::Normalize(l2) >= 0.0 &&
             vtkMath::Dot(l1,l2) < CosEdgeAngle)
        {
          node_type->SetValue(id_node, VTK_FIXED_VERTEX);
        }
      }//if along edge
    }//if edge vertex
  }
  cout<<"m_PotentialSnapPoints.size()="<<m_PotentialSnapPoints.size()<<endl;
  cout<<"=== UpdateNodeType END ==="<<endl;
  return(0);
}

///@@@  TODO: Optimize
char SurfaceOperation::getNodeType(vtkIdType id_node)
{
  //initialize default value
  char type=VTK_SIMPLE_VERTEX;
  
  //loop through edges around id_node
  
  QVector <vtkIdType> edges;
  
  double CosEdgeAngle = cos((double) vtkMath::RadiansFromDegrees(this->EdgeAngle));
  
  foreach(int i_node2, n2n[_nodes[id_node]])
  {
    vtkIdType id_node2 = nodes[i_node2];
    //-----------------------
    //determine edge type
    char edge = getEdgeType(id_node2,id_node);
    
    //-----------------------
    //determine node type pre-processing (count nb of complex edges if the node is complex, otherwise, just count the nb of edges)
    if ( edge && type == VTK_SIMPLE_VERTEX )
    {
      edges.clear();
      edges.push_back(id_node2);
      type = edge;
    }
    else if ( (edge && type == VTK_BOUNDARY_EDGE_VERTEX) ||
              (edge && type == VTK_FEATURE_EDGE_VERTEX) ||
              (!edge && type == VTK_SIMPLE_VERTEX ) )
    {
      edges.push_back(id_node2);
      if ( type && edge == VTK_BOUNDARY_EDGE_VERTEX )
      {
        type = VTK_BOUNDARY_EDGE_VERTEX;//VTK_BOUNDARY_EDGE_VERTEX has priority over VTK_FEATURE_EDGE_VERTEX
      }
    }
  }
  //-----------------------
  //determine node type post-processing
  if ( type == VTK_FEATURE_EDGE_VERTEX || type == VTK_BOUNDARY_EDGE_VERTEX )
  { //see how many edges; if two, what the angle is
    
    if ( !this->BoundarySmoothing && 
         type == VTK_BOUNDARY_EDGE_VERTEX )
    {
      type = VTK_FIXED_VERTEX;
    }
    else if ( edges.size() != 2 )
    {
      type = VTK_FIXED_VERTEX;
    }
    else //check angle between edges
    {
      double x1[3], x2[3], x3[3], l1[3], l2[3];
      grid->GetPoint(edges[0],x1);
      grid->GetPoint(id_node,x2);
      grid->GetPoint(edges[1],x3);
      for (int k=0; k<3; k++)
      {
        l1[k] = x2[k] - x1[k];
        l2[k] = x3[k] - x2[k];
      }
      if ( vtkMath::Normalize(l1) >= 0.0 &&
           vtkMath::Normalize(l2) >= 0.0 &&
           vtkMath::Dot(l1,l2) < CosEdgeAngle)
      {
        type = VTK_FIXED_VERTEX;
      }
    }//if along edge
  }//if edge vertex
  
  return(type);
}

int SurfaceOperation::getEdgeCells(vtkIdType id_node1, vtkIdType id_node2,QVector <vtkIdType> &EdgeCells)
{
  QSet<vtkIdType> S1;
  foreach(int i, n2c[_nodes[id_node1]]) S1.insert(cells[i]);
  
  QSet<vtkIdType> S2;
  foreach(int i, n2c[_nodes[id_node2]]) S2.insert(cells[i]);
  
  S2.intersect(S1);
  EdgeCells = Set2Vector(S2,false);
  return EdgeCells.size();
}

int SurfaceOperation::getEdgeCells(vtkIdType id_node1, vtkIdType id_node2,QSet <vtkIdType> &EdgeCells)
{
  QSet<vtkIdType> S1;
  foreach(int i, n2c[_nodes[id_node1]]) S1.insert(cells[i]);
  
  QSet<vtkIdType> S2;
  foreach(int i, n2c[_nodes[id_node2]]) S2.insert(cells[i]);
  
  EdgeCells = S2.intersect(S1);
  return EdgeCells.size();
}

char SurfaceOperation::getEdgeType(vtkIdType a_node1, vtkIdType a_node2)
{
  double CosFeatureAngle = cos((double) vtkMath::RadiansFromDegrees(this->FeatureAngle));
  
  //compute number of cells around edge [a_node,p2] and put them into neighbour_cells
  QVector <vtkIdType> neighbour_cells;
  int numNei = getEdgeCells(a_node1,a_node2,neighbour_cells) - 1;
  
  //set default value
  char edge = VTK_SIMPLE_VERTEX;
  
  if ( numNei == 0 )
  {
    edge = VTK_BOUNDARY_EDGE_VERTEX;
  }
  else if ( numNei >= 2 )
  {
    edge = VTK_FEATURE_EDGE_VERTEX;
  }
  else if ( numNei == 1 )
  {
    //check angle between cell1 and cell2 against FeatureAngle
    if ( this->FeatureEdgeSmoothing && CosAngle(grid,neighbour_cells[0],neighbour_cells[1]) <= CosFeatureAngle )
    {
      edge = VTK_FEATURE_EDGE_VERTEX;
    }
    //check the boundary codes
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    if ( cell_code->GetValue(neighbour_cells[0]) !=  cell_code->GetValue(neighbour_cells[1]) )
    {
      edge = VTK_BOUNDARY_EDGE_VERTEX;
    }
  }
  
  return(edge);
}

QSet <int> SurfaceOperation::getBCset(vtkIdType id_node)
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  QSet <int> bc;
  foreach(int i_cell, n2c[_nodes[id_node]])
  {
    vtkIdType id_cell = cells[i_cell];
    bc.insert(cell_code->GetValue(id_cell));
  }
  return(bc);
}

VertexMeshDensity SurfaceOperation::getVMD(vtkIdType id_node)
{
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  
  VertexMeshDensity VMD;
  VMD.type = node_type->GetValue(id_node);
  VMD.density = 0;
  VMD.CurrentNode = id_node;
  
  foreach(int i_cell, n2c[_nodes[id_node]])
  {
    vtkIdType id_cell = cells[i_cell];
    VMD.BCmap[cell_code->GetValue(id_cell)] = 2;
  }
  return(VMD);
}

void SurfaceOperation::setSource(vtkUnstructuredGrid *a_ProjectionSurface)
{
  if(m_CellLocator) {
    cout<<"WARNING: Deleting previous m_CellLocator!"<<endl;
    m_CellLocator->Delete();
    m_CellLocator=NULL;
  }
  if(m_ProjectionSurface) {
    cout<<"WARNING: Deleting previous m_ProjectionSurface!"<<endl;
    m_ProjectionSurface->Delete();
    m_ProjectionSurface=NULL;
  }
  
  m_ProjectionSurface=vtkUnstructuredGrid::New();
  makeCopy(a_ProjectionSurface,m_ProjectionSurface);
  
  m_CellLocator=vtkCellLocator::New();
  m_CellLocator->SetDataSet(a_ProjectionSurface);
  m_CellLocator->BuildLocator();
//   m_CellLocator->CacheCellBoundsOn();
  cout<<"m_CellLocator->GetNumberOfBuckets()="<<m_CellLocator->GetNumberOfBuckets()<<endl;
  cout<<"m_CellLocator->GetNumberOfCellsPerBucket()="<<m_CellLocator->GetNumberOfCellsPerBucket()<<endl;
  cout<<"m_CellLocator->GetCacheCellBounds()="<<m_CellLocator->GetCacheCellBounds()<<endl;
  
  cout<<"ORIGINAL: m_CellLocator="<<m_CellLocator<<endl;
  cout<<"ORIGINAL: m_ProjectionSurface="<<m_ProjectionSurface<<endl;
}

void SurfaceOperation::set_CellLocator_and_ProjectionSurface(vtkCellLocator *a_CellLocator, vtkUnstructuredGrid *a_ProjectionSurface)
{
  m_CellLocator = vtkCellLocator::SafeDownCast(a_CellLocator);
  m_ProjectionSurface = vtkUnstructuredGrid::SafeDownCast(a_ProjectionSurface);
  
  cout<<"===set_CellLocator_and_ProjectionSurface==="<<endl;
  cout_grid(cout,m_ProjectionSurface);
  
  cout<<"COPY: m_CellLocator="<<m_CellLocator<<endl;
  cout<<"COPY: m_ProjectionSurface="<<m_ProjectionSurface<<endl;
}

vec3_t SurfaceOperation::project(vec3_t OM)
{
  vec3_t OP;
  if(m_CellLocator==NULL) {
    cout<<"FATAL ERROR: No source surface has been defined."<<endl; EG_BUG;
  }
  else {
    vtkIdType cellId;
    int subId;
    double dist2;
    m_CellLocator->FindClosestPoint(OM.data(),OP.data(),cellId,subId,dist2);
  }
  
//   OM=OA+AP+PM;
/*  vec3_t OA(0,0,OM[2]);
  vec3_t AM=OM-OA;
  vec3_t r=AM;
  r.normalise();
  vec3_t AP=0.1*r;
  OP=OA+AP;*/
  
  return(OP);
}

void SurfaceOperation::delete_CellLocator_and_ProjectionSurface()
{
  if(m_CellLocator) {
    cout<<"WARNING: Deleting m_CellLocator!"<<endl;
    m_CellLocator->Delete();
    m_CellLocator=NULL;
  }
  if(m_ProjectionSurface) {
    cout<<"WARNING: Deleting m_ProjectionSurface!"<<endl;
    m_ProjectionSurface->Delete();
    m_ProjectionSurface=NULL;
  }
}

//////////////////////////////////////////////
double SurfaceOperation::CurrentVertexAvgDist(vtkIdType id_node)
{
  double total_dist=0;
  double avg_dist=0;
  int N = n2n[_nodes[id_node]].size();
  vec3_t C;
  grid->GetPoint(id_node, C.data());
  foreach(int i_node_neighbour, n2n[_nodes[id_node]])
  {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    vec3_t M;
    grid->GetPoint(id_node_neighbour, M.data());
    total_dist+=(M-C).abs();
  }
  avg_dist=total_dist/(double)N;
  return(avg_dist);
}

double SurfaceOperation::CurrentMeshDensity(vtkIdType id_node)
{
  return 1./CurrentVertexAvgDist(id_node);
}

double SurfaceOperation::DesiredVertexAvgDist(vtkIdType id_node)
{
  double total_dist=0;
  double avg_dist=0;
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  int N = n2n[_nodes[id_node]].size();
  foreach(int i_node_neighbour, n2n[_nodes[id_node]])
  {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    total_dist+=1./node_meshdensity_desired->GetValue(id_node_neighbour);
  }
  avg_dist=total_dist/(double)N;
  return(avg_dist);
}

double SurfaceOperation::DesiredMeshDensity(vtkIdType id_node)
{
  double total_density=0;
  double avg_density=0;
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  int N = n2n[_nodes[id_node]].size();
  foreach(int i_node_neighbour, n2n[_nodes[id_node]])
  {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    total_density+=node_meshdensity_desired->GetValue(id_node_neighbour);
  }
  avg_density=total_density/(double)N;
  return(avg_density);
}

///////////////////////////////////////////

//---------------------------------------------------
//Utility functions used in Roland's formulas
//Should be renamed to be more explicit
//Some could be moved into geometrytools
//Some are pretty useless

///perimeter
double SurfaceOperation::Um(vtkIdType D) {
  double ret=0;
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(D, N_pts, pts);
  for(int i=0;i<N_pts;i++)
  {
    vec3_t A,B;
    grid->GetPoints()->GetPoint(pts[i], A.data());
    grid->GetPoints()->GetPoint(pts[(i+1)%N_pts], B.data());
    ret+=(B-A).abs();
  }
  return(ret);
}

/// area of the circumscribed circle of the triangle
double SurfaceOperation::A_U(vtkIdType D) {
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(D, N_pts, pts);
  vec3_t A,B,C;
  grid->GetPoints()->GetPoint(pts[0], A.data());
  grid->GetPoints()->GetPoint(pts[1], B.data());
  grid->GetPoints()->GetPoint(pts[2], C.data());
  double a=(C-B).abs();
  double alpha=angle((B-A),(C-A));
  double R=a/(2*sin(alpha));
  return(M_PI*R*R);
}

/// triangle area
double SurfaceOperation::A_D(vtkIdType D) {
  return(cellVA(grid,D));
}

/// triangle neighbours
double SurfaceOperation::DN(int i,vtkIdType D) {
  return(c2c[D][i]);
}

/// number of edges
double SurfaceOperation::nk(vtkIdType P) {
  return(n2n[P].size());
}

double SurfaceOperation::G_k(vtkIdType node) {
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  return(1.0/node_meshdensity_desired->GetValue(node));
}

/// triangle nodes
double SurfaceOperation::DK(int i,vtkIdType D) {
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(D, N_pts, pts);
  return(pts[i]);
}

vtkIdType SurfaceOperation::KK(int i,vtkIdType j,vtkIdType K) {//i=1 or 2, j=node2, K=node1
  if(i==1) return(K);
  else return(j);
}

double SurfaceOperation::L_k(vtkIdType j,vtkIdType K)// node1 K, node2 j
{
  vec3_t A;
  vec3_t B;
  grid->GetPoints()->GetPoint(K, A.data());
  grid->GetPoints()->GetPoint(j, B.data());
  return((B-A).abs());
}

double SurfaceOperation::Q_L(vtkIdType D)
{
      // Um(D)/sum(G_k(DK(i,D)),i,1,3)
  double denom_sum=0;
  for(int i=0;i<3;i++)
  {
    denom_sum += G_k(DK(i,D));
  }
      /*if(DebugLevel>0) cout<<"D="<<D<<" Um(D)="<<Um(D)<<" denom_sum="<<denom_sum<<endl;*/
  return(Um(D)/denom_sum);
}

double SurfaceOperation::Q_L1(vtkIdType P)
{
      // [2*sum(L_k(i~),i,1,nk(P))]/[sum(G_k(KK(1,i~))+G_k(KK(2,i~)),i,1,nk(P))]
  double num_sum=0;
  double denom_sum=0;
  foreach(vtkIdType j,n2n[P])
  {
    num_sum += 2*L_k(j,P);
    denom_sum += G_k(KK(1,j,P))+G_k(KK(2,j,P));
  }
  return(num_sum/denom_sum);
}

double SurfaceOperation::Q_L2(vtkIdType P)
{
      // min([2*L_k(i~)]/[G_k(KK(1,i~))+G_k(KK(2,i~))])
  QVector <double> V;
  double num,denom;
  foreach(vtkIdType j,n2n[P])
  {
    num = 2*L_k(j,P);
    denom = G_k(KK(1,j,P))+G_k(KK(2,j,P));
    V.push_back(num/denom);
  }
  qSort(V.begin(),V.end());
  return(V[0]);
}

double SurfaceOperation::T_min(int w)
{
  // sum([A_U(i)]/[A_D(i)^w]*[G_k(i)^(2*(w-1))],i,1,Nd)
  int N_cells=grid->GetNumberOfCells();
  double T=0;
  for(int i=0;i<N_cells;i++)
  {
    T += A_U(i)/pow(A_D(i),w)*pow(G_k(i),2*(w-1));
  }
  return(T);
}

//---------------------------------------------------

vtkIdType SurfaceOperation::getClosestNode(vtkIdType id_node)
{
  vec3_t C;
  grid->GetPoint(id_node,C.data());
  vtkIdType id_minlen=-1;
  double minlen=-1;
  foreach(vtkIdType neighbour,n2n[id_node])
  {
    vec3_t M;
    grid->GetPoint(neighbour,M.data());
    double len=(M-C).abs();
    if(minlen<0 or len<minlen)
    {
      minlen=len;
      id_minlen=neighbour;
    }
  }
  return(id_minlen);
}

vtkIdType SurfaceOperation::getFarthestNode(vtkIdType id_node)
{
  vec3_t C;
  grid->GetPoint(id_node,C.data());
  vtkIdType id_maxlen=-1;
  double maxlen=-1;
  foreach(vtkIdType neighbour,n2n[id_node])
  {
    vec3_t M;
    grid->GetPoint(neighbour,M.data());
    double len=(M-C).abs();
    if(maxlen<0 or len>maxlen)
    {
      maxlen=len;
      id_maxlen=neighbour;
    }
  }
  return(id_maxlen);
}

int SurfaceOperation::NumberOfCommonPoints(vtkIdType node1, vtkIdType node2, bool& IsTetra)
{
//   QVector< QSet< int > > 	n2n
  QSet <int> node1_neighbours=n2n[node1];
  QSet <int> node2_neighbours=n2n[node2];
  QSet <int> intersection=node1_neighbours.intersect(node2_neighbours);
  int N=intersection.size();
  IsTetra=false;
  if(N==2)
  {
    QSet<int>::const_iterator p1=intersection.begin();
    QSet<int>::const_iterator p2=p1+1;
    vtkIdType intersection1=_nodes[*p1];
    vtkIdType intersection2=_nodes[*p2];
    if(n2n[intersection1].contains(intersection2))//if there's an edge between intersection1 and intersection2
    {
      //check if (node1,intersection1,intersection2) and (node2,intersection1,intersection2) are defined as cells!
  //     QVector< QSet< int > > 	n2c
      QSet< int > S1=n2c[intersection1];
      QSet< int > S2=n2c[intersection2];
      QSet< int > Si=S1.intersect(S2);
      int counter=0;
      foreach(vtkIdType C,Si){
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(C, N_pts, pts);
        for(int i=0;i<N_pts;i++)
        {
          if(pts[i]==node1 || pts[i]==node2) counter++;
        }
      }
      if(counter>=2) IsTetra=true;
    }
  }
  return(N);
}

QVector <vtkIdType> SurfaceOperation::getPotentialSnapPoints(vtkIdType id_node)
{
  if(id_node<0 || id_node>=m_PotentialSnapPoints.size()) EG_BUG;
  return m_PotentialSnapPoints[id_node];
}

// DEFINITIONS:
// Normal cell: nothing has changed
// Dead cell: the cell does not exist anymore
// Mutated cell: the cell's form has changed
// Mutilated cell: the cell has less points than before

vtkIdType SurfaceOperation::FindSnapPoint(vtkIdType DeadNode,QSet <vtkIdType> & DeadCells,QSet <vtkIdType> & MutatedCells,QSet <vtkIdType> & MutilatedCells, int& N_newpoints, int& N_newcells)
{
  ///@@@  TODO: Organize cases and make sure all are considered if possible.
  getAllSurfaceCells(cells,grid);
  getNodesFromCells(cells, nodes, grid);
  setGrid(grid);
  setCells(cells);
  
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  if(node_type->GetValue(DeadNode)==VTK_FIXED_VERTEX )
  {
    cout<<"Sorry, unable to remove fixed vertex."<<endl;
    return(-1);
  }
  
  //grid info
  int N_points=grid->GetNumberOfPoints();
  int N_cells=grid->GetNumberOfCells();
  N_newpoints=-1;
  N_newcells=0;
  
  vtkIdType SnapPoint=-1;
  
  foreach(vtkIdType PSP, n2n[DeadNode])
  {
    bool IsValidSnapPoint=true;
    
    if(DebugLevel>10) cout<<"====>PSP="<<PSP<<endl;
    bool IsTetra=true;
    if(NumberOfCommonPoints(DeadNode,PSP,IsTetra)>2)//common point check
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    if(IsTetra)//tetra check
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    //count number of points and cells to remove + analyse cell transformations
    N_newpoints=-1;
    N_newcells=0;
    DeadCells.clear();
    MutatedCells.clear();
    MutilatedCells.clear();
    foreach(vtkIdType C, n2c[DeadNode])//loop through potentially dead cells
    {
      //get points around cell
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(C, N_pts, pts);
      
      bool ContainsSnapPoint=false;
      bool invincible=false;
      for(int i=0;i<N_pts;i++)
      {
        if(DebugLevel>10) cout<<"pts["<<i<<"]="<<pts[i]<<" and PSP="<<PSP<<endl;
        if(pts[i]==PSP) {ContainsSnapPoint=true;}
        if(pts[i]!=DeadNode && pts[i]!=PSP &&  n2c[pts[i]].size()<=1) invincible=true;
      }
      if(ContainsSnapPoint)
      {
        if(N_pts==3)//dead cell
        {
          if(invincible)//Check that empty lines aren't left behind when a cell is killed
          {
            IsValidSnapPoint=false;
          }
          else
          {
            DeadCells.insert(C);
            N_newcells-=1;
          }
        }
        else
        {
          cout<<"RED ALERT: Xenomorph detected!"<<endl;
          EG_BUG;
        }
      }
      else
      {
        vtkIdType src_N_pts, *src_pts;
        grid->GetCellPoints(C, src_N_pts, src_pts);
        
        if(src_N_pts!=3)
        {
          cout<<"RED ALERT: Xenomorph detected!"<<endl;
          EG_BUG;
        }
        
        vtkIdType OldTriangle[3];
        vtkIdType NewTriangle[3];
        
        for(int i=0;i<src_N_pts;i++)
        {
          OldTriangle[i]=src_pts[i];
          NewTriangle[i]=( (src_pts[i]==DeadNode) ? PSP : src_pts[i] );
        }
        vec3_t Old_N= triNormal(grid, OldTriangle[0], OldTriangle[1], OldTriangle[2]);
        vec3_t New_N= triNormal(grid, NewTriangle[0], NewTriangle[1], NewTriangle[2]);
        
        if(Old_N*New_N<Old_N*Old_N*1./100.)//area + inversion check
        {
          if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
          IsValidSnapPoint=false;
        }
        
        //mutated cell
        MutatedCells.insert(C);
        if(DebugLevel>10) cout<<"cell "<<C<<" has been infected!"<<endl;
      }
    }
    
    if(N_cells+N_newcells<=0)//survivor check
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    if(node_type->GetValue(DeadNode)==VTK_BOUNDARY_EDGE_VERTEX && node_type->GetValue(PSP)==VTK_SIMPLE_VERTEX)
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    if(node_type->GetValue(DeadNode)==VTK_BOUNDARY_EDGE_VERTEX)
    {
      QVector <vtkIdType> Peons = getPotentialSnapPoints(DeadNode);
      if(!Peons.contains(PSP))
      {
        if(DebugLevel>0) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
        IsValidSnapPoint=false;
      }
    }
    
    if(node_type->GetValue(DeadNode)==VTK_FEATURE_EDGE_VERTEX && node_type->GetValue(PSP)==VTK_SIMPLE_VERTEX)
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    if(node_type->GetValue(DeadNode)==VTK_FEATURE_EDGE_VERTEX)
    {
      QVector <vtkIdType> Peons = getPotentialSnapPoints(DeadNode);
      if(!Peons.contains(PSP))
      {
        if(DebugLevel>0) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
        IsValidSnapPoint=false;
      }
    }
    
    if(IsValidSnapPoint) {SnapPoint=PSP; break;}
  }//end of loop through potential SnapPoints
  
  if(DebugLevel>10)
  {
    cout<<"AT FINDSNAPPOINT EXIT"<<endl;
    cout<<"N_points="<<N_points<<endl;
    cout<<"N_cells="<<N_cells<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
  }
  return(SnapPoint);
}
//End of FindSnapPoint

bool SurfaceOperation::DeletePoint(vtkIdType DeadNode, int& N_newpoints, int& N_newcells)
{
  QSet <vtkIdType> DeadNodes;
  DeadNodes.insert(DeadNode);
  bool ret = DeleteSetOfPoints(DeadNodes, N_newpoints, N_newcells);
  return(ret);
}
//End of DeletePoint

bool SurfaceOperation::DeleteSetOfPoints(QSet <vtkIdType> DeadNodes, int& N_newpoints, int& N_newcells)
{
  QVector <vtkIdType> DeadNode_vector=Set2Vector(DeadNodes,false);
  
  getAllSurfaceCells(cells,grid);
  UpdateNodeType();
  
  //src grid info
  int N_points = grid->GetNumberOfPoints();
  int N_cells = grid->GetNumberOfCells();
  
  QSet <vtkIdType> DeadCells;
  QSet <vtkIdType> MutatedCells;
  QSet <vtkIdType> MutilatedCells;
  QVector <vtkIdType> SnapPoint(DeadNode_vector.size());
  
  //counter init
  N_newpoints=0;
  N_newcells=0;
  
  for(int i=0;i<DeadNode_vector.size();i++)
  {
    if(DeadNode_vector[i]<0 || DeadNode_vector[i]>=N_points)
    {
      cout<<"Warning: Point out of range: DeadNode_vector[i]="<<DeadNode_vector[i]<<" N_points="<<N_points<<endl;
      return(false);
    }
    
    //local values
    int l_N_newpoints;
    int l_N_newcells;
    QSet <vtkIdType> l_DeadCells;
    QSet <vtkIdType> l_MutatedCells;
    QSet <vtkIdType> l_MutilatedCells;
    
    SnapPoint[i]=FindSnapPoint(DeadNode_vector[i], l_DeadCells, l_MutatedCells, l_MutilatedCells, l_N_newpoints, l_N_newcells);
    //global values
    N_newpoints+=l_N_newpoints;
    N_newcells+=l_N_newcells;
    DeadCells.unite(l_DeadCells);//DeadCells unite! Kill the living! :D
    MutatedCells.unite(l_MutatedCells);
    MutilatedCells.unite(l_MutilatedCells);
    
    if(DebugLevel>0) cout<<"===>DeadNode_vector[i]="<<DeadNode_vector[i]<<" moving to SNAPPOINT="<<SnapPoint[i]<<" DebugLevel="<<DebugLevel<<endl;
    if(SnapPoint[i]<0) {cout<<"Sorry no possible SnapPoint found."<<endl; return(false);}
    
  }
  
  //allocate
  EG_VTKSP(vtkUnstructuredGrid,dst);
  allocateGrid(dst,N_cells+N_newcells,N_points+N_newpoints);
  
  //vector used to redefine the new point IDs
  QVector <vtkIdType> OffSet(N_points);
  
  //copy undead points
  vtkIdType dst_id_node=0;
  for (vtkIdType src_id_node = 0; src_id_node < N_points; src_id_node++) {//loop through src points
    if(!DeadNode_vector.contains(src_id_node))//if the node isn't dead, copy it
    {
      vec3_t x;
      grid->GetPoints()->GetPoint(src_id_node, x.data());
      dst->GetPoints()->SetPoint(dst_id_node, x.data());
      copyNodeData(grid, src_id_node, dst, dst_id_node);
      OffSet[src_id_node]=src_id_node-dst_id_node;
      dst_id_node++;
    }
    else
    {
      if(DebugLevel>0) cout<<"src_id_node="<<src_id_node<<" dst_id_node="<<dst_id_node<<endl;
    }
  };
  //Copy undead cells
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {//loop through src cells
    if(!DeadCells.contains(id_cell))//if the cell isn't dead
    {
      vtkIdType src_N_pts, *src_pts;
      vtkIdType dst_N_pts, *dst_pts;
      grid->GetCellPoints(id_cell, src_N_pts, src_pts);
      
      vtkIdType type_cell = grid->GetCellType(id_cell);
//       src->GetCellPoints(id_cell, dst_N_pts, dst_pts);
      dst_N_pts=src_N_pts;
      dst_pts=new vtkIdType[dst_N_pts];
      if(MutatedCells.contains(id_cell))//mutated cell
      {
        for(int i=0;i<src_N_pts;i++)
        {
          int DeadIndex = DeadNode_vector.indexOf(src_pts[i]);
          if(DeadIndex!=-1) {
            dst_pts[i]=SnapPoint[DeadIndex]-OffSet[SnapPoint[DeadIndex]];
          }
          else dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
      }
      else if(MutilatedCells.contains(id_cell))//mutilated cell (ex: square becoming triangle) (WARNING: Not fully functional yet)
      {
        cout<<"FATAL ERROR: Quads not supported yet."<<endl;EG_BUG;
        
        if(type_cell==VTK_QUAD) {
          type_cell=VTK_TRIANGLE;
          dst_N_pts-=1;
        }
        else {cout<<"FATAL ERROR: Unknown mutilated cell detected! It is not a quad! Potential xenomorph infestation!"<<endl;EG_BUG;}
        //merge points
        for(int i=0;i<src_N_pts;i++)
        {
/*          if(src_pts[i]==SnapPoint) { dst_pts[j]=SnapPoint-OffSet[SnapPoint];j++; }//SnapPoint
          else if(src_pts[i]!=DeadNode_vector[i]) { dst_pts[j]=src_pts[i]-OffSet[src_pts[i]];j++; }//pre-snap/dead + post-snap/dead*/
          //do nothing in case of DeadNode_vector[i]
        }
      }
      else//normal cell
      {
        if(DebugLevel>10) cout<<"processing normal cell "<<id_cell<<endl;
        for(int i=0;i<src_N_pts;i++)
        {
          dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
      }
      //copy the cell
      vtkIdType id_new_cell = dst->InsertNextCell(type_cell, dst_N_pts, dst_pts);
      copyCellData(grid, id_cell, dst, id_new_cell);
      delete dst_pts;
    }
  };
  
  makeCopy(dst, grid);
  return(true);
}
//End of DeleteSetOfPoints

bool SurfaceOperation::FlippedCells(vtkIdType id_node, vec3_t P)
{
  vec3_t x0_old, x0_new;
  grid->GetPoint(id_node, x0_old.data());
  x0_new=P;
  
  foreach(int i_cell, n2c[_nodes[id_node]])
  {
    vtkIdType id_cell = cells[i_cell];
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    int i;
    for(i=0;i<N_pts;i++)
    {
      if(pts[i]==id_node) break;
    }
    vec3_t x2, x3;
    grid->GetPoint(pts[(i+1)%N_pts], x2.data());
    grid->GetPoint(pts[(i+2)%N_pts], x3.data());
    vec3_t v2_old=x2-x0_old;
    vec3_t v3_old=x3-x0_old;
    
    //top point
    vec3_t S=v2_old.cross(v3_old);
    double V_old=tetraVol(x0_old, S, x2, x3, true);
    double V_new=tetraVol(x0_new, S, x2, x3, true);
    double prod=V_old*V_new;
    if( prod<0 ) {
      return(true);
    }
  }
  return(false);
}