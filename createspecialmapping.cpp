//
// C++ Implementation: createspecialmapping
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "createspecialmapping.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTriangleFilter.h"
#include <QString>
#include <QTextStream>

#include <iostream>
using namespace std;

#define VTK_SIMPLE_VERTEX 0
#define VTK_FIXED_VERTEX 1
#define VTK_FEATURE_EDGE_VERTEX 2
#define VTK_BOUNDARY_EDGE_VERTEX 3

const char* VertexType(char T)
{
  if(T==VTK_SIMPLE_VERTEX) return("VTK_SIMPLE_VERTEX");
  if(T==VTK_FIXED_VERTEX) return("VTK_FIXED_VERTEX");
  if(T==VTK_FEATURE_EDGE_VERTEX) return("VTK_FEATURE_EDGE_VERTEX");
  if(T==VTK_BOUNDARY_EDGE_VERTEX) return("VTK_BOUNDARY_EDGE_VERTEX");
  else return("Unknown vertex type");
}

// Special structure for marking vertices
typedef struct _vtkMeshVertex 
{
  char      type;
  vtkIdList *edges; // connected edges (list of connected point ids)
} vtkMeshVertex, *vtkMeshVertexPtr;

CreateSpecialMapping::CreateSpecialMapping()
{
}

int CreateSpecialMapping::Process(vtkPolyData* input)
{
  cout<<"=================="<<endl;
  
  vtkPolyData *source = 0;
  
  vtkIdType numPts, numCells, i, numPolys, numStrips;
  int j, k;
  vtkIdType npts = 0;
  vtkIdType *pts = 0;
  vtkIdType p1, p2;
  double x[3], y[3], deltaX[3], xNew[3], conv, maxDist, dist, factor;
  double x1[3], x2[3], x3[3], l1[3], l2[3];
  double CosFeatureAngle; //Cosine of angle between adjacent polys
  double CosEdgeAngle; // Cosine of angle between adjacent edges
  double closestPt[3], dist2, *w = NULL;
  int iterationNumber, abortExecute;
  vtkIdType numSimple=0, numBEdges=0, numFixed=0, numFEdges=0;
  vtkPolyData *inMesh, *Mesh;
  vtkPoints *inPts;
  vtkTriangleFilter *toTris=NULL;
  vtkCellArray *inVerts, *inLines, *inPolys, *inStrips;
  vtkPoints *newPts;
  vtkMeshVertexPtr Verts;
  vtkCellLocator *cellLocator=NULL;
  
  // Check input
  //
  numPts=input->GetNumberOfPoints();
  numCells=input->GetNumberOfCells();
  if (numPts < 1 || numCells < 1)
  {
    cout<<"No data to smooth!"<<endl;
    return 1;
  }
  
  CosFeatureAngle = 
    cos((double) vtkMath::DegreesToRadians() * this->FeatureAngle);
  CosEdgeAngle = cos((double) vtkMath::DegreesToRadians() * this->EdgeAngle);
  
  cout<<"Smoothing " << numPts << " vertices, " << numCells 
                << " cells with:\n"
                << "\tConvergence= " << this->Convergence << "\n"
                << "\tIterations= " << this->NumberOfIterations << "\n"
                << "\tRelaxation Factor= " << this->RelaxationFactor << "\n"
                << "\tEdge Angle= " << this->EdgeAngle << "\n"
                << "\tBoundary Smoothing " << (this->BoundarySmoothing ? "On\n" : "Off\n")
                << "\tFeature Edge Smoothing " << (this->FeatureEdgeSmoothing ? "On\n" : "Off\n")
                << "\tError Scalars " << (this->GenerateErrorScalars ? "On\n" : "Off\n")
                << "\tError Vectors " << (this->GenerateErrorVectors ? "On\n" : "Off\n")<<endl;
  
  // Peform topological analysis. What we're gonna do is build a connectivity
  // array of connected vertices. The outcome will be one of three
  // classifications for a vertex: VTK_SIMPLE_VERTEX, VTK_FIXED_VERTEX. or
  // VTK_EDGE_VERTEX. Simple vertices are smoothed using all connected 
  // vertices. FIXED vertices are never smoothed. Edge vertices are smoothed
  // using a subset of the attached vertices.
  //
  cout<<"Analyzing topology..."<<endl;
  cout<<"0:numPts="<<numPts<<endl;
  Verts = new vtkMeshVertex[numPts];
  for (i=0; i<numPts; i++)
  {
    cout<<"0:VTK_SIMPLE_VERTEX"<<endl;
    Verts[i].type = VTK_SIMPLE_VERTEX; //can smooth
    Verts[i].edges = NULL;
  }
  
  inPts = input->GetPoints();
  
  // check vertices first. Vertices are never smoothed_--------------
  for (inVerts=input->GetVerts(), inVerts->InitTraversal(); 
       inVerts->GetNextCell(npts,pts); )
  {
    for (j=0; j<npts; j++)
    {
      cout<<"pts[j]="<<pts[j]<<"->vertices:VTK_FIXED_VERTEX"<<endl;
      Verts[pts[j]].type = VTK_FIXED_VERTEX;
    }
  }
  
  // now check lines. Only manifold lines can be smoothed------------
  for (inLines=input->GetLines(), inLines->InitTraversal(); 
       inLines->GetNextCell(npts,pts); )
  {
    for (j=0; j<npts; j++)
    {
      cout<<"pts[j]="<<pts[j]<<"->lines"<<endl;
      if ( Verts[pts[j]].type == VTK_SIMPLE_VERTEX )
      {
        if ( j == (npts-1) ) //end-of-line marked FIXED
        {
          cout<<"pts[j]="<<pts[j]<<"2:VTK_FIXED_VERTEX"<<endl;
          Verts[pts[j]].type = VTK_FIXED_VERTEX;
        }
        else if ( j == 0 ) //beginning-of-line marked FIXED
        {
          cout<<"pts[j]="<<pts[j]<<"3:VTK_FIXED_VERTEX"<<endl;
          Verts[pts[0]].type = VTK_FIXED_VERTEX;
          inPts->GetPoint(pts[0],x2);
          inPts->GetPoint(pts[1],x3);
        }
        else //is edge vertex (unless already edge vertex!)
        {
          cout<<"pts[j]="<<pts[j]<<"4:VTK_FEATURE_EDGE_VERTEX"<<endl;
          Verts[pts[j]].type = VTK_FEATURE_EDGE_VERTEX;
          Verts[pts[j]].edges = vtkIdList::New();
          Verts[pts[j]].edges->SetNumberOfIds(2);
          Verts[pts[j]].edges->SetId(0,pts[j-1]);
          Verts[pts[j]].edges->SetId(1,pts[j+1]);
        }
      } //if simple vertex
      
      else if ( Verts[pts[j]].type == VTK_FEATURE_EDGE_VERTEX )
      { //multiply connected, becomes fixed!
        cout<<"pts[j]="<<pts[j]<<"5:VTK_FIXED_VERTEX"<<endl;
        Verts[pts[j]].type = VTK_FIXED_VERTEX;
        Verts[pts[j]].edges->Delete();
        Verts[pts[j]].edges = NULL;
      }
      
    } //for all points in this line
  } //for all lines
  
  // now polygons and triangle strips-------------------------------
  inPolys=input->GetPolys();
  numPolys = inPolys->GetNumberOfCells();
  inStrips=input->GetStrips();
  numStrips = inStrips->GetNumberOfCells();
  
  cout<<"numPolys="<<numPolys<<endl;
  cout<<"numStrips="<<numStrips<<endl;
  
  
  if ( numPolys > 0 || numStrips > 0 )
  { //build cell structure
    vtkCellArray *polys;
    vtkIdType cellId;
    int numNei, nei, edge;
    vtkIdType numNeiPts;
    vtkIdType *neiPts;
    double normal[3], neiNormal[3];
    vtkIdList *neighbors;
    
    neighbors = vtkIdList::New();
    neighbors->Allocate(VTK_CELL_SIZE);
    
    inMesh = vtkPolyData::New();
    inMesh->SetPoints(inPts);
    inMesh->SetPolys(inPolys);
    Mesh = inMesh;
    
    if ( (numStrips = inStrips->GetNumberOfCells()) > 0 )
    { // convert data to triangles
      inMesh->SetStrips(inStrips);
      toTris = vtkTriangleFilter::New();
      toTris->SetInput(inMesh);
      toTris->Update();
      Mesh = toTris->GetOutput();
    }
    
    Mesh->BuildLinks(); //to do neighborhood searching
    polys = Mesh->GetPolys();
    
    for (cellId=0, polys->InitTraversal(); polys->GetNextCell(npts,pts); 
         cellId++)
    {
      cout<<"->cellId="<<cellId<<endl;
      for (i=0; i < npts; i++) 
      {
        cout<<"-->i="<<i<<endl;
        p1 = pts[i];
        p2 = pts[(i+1)%npts];
        
        if ( Verts[p1].edges == NULL )
        {
          Verts[p1].edges = vtkIdList::New();
          Verts[p1].edges->Allocate(16,6);
        }
        if ( Verts[p2].edges == NULL )
        {
          Verts[p2].edges = vtkIdList::New();
          Verts[p2].edges->Allocate(16,6);
        }
        
        Mesh->GetCellEdgeNeighbors(cellId,p1,p2,neighbors);
        numNei = neighbors->GetNumberOfIds();
        cout<<"-->numNei="<<numNei<<endl;
        
        edge = VTK_SIMPLE_VERTEX;
        if ( numNei == 0 )
        {
          edge = VTK_BOUNDARY_EDGE_VERTEX;
        }
        
        else if ( numNei >= 2 )
        {
          // check to make sure that this edge hasn't been marked already
          for (j=0; j < numNei; j++)
          {
            if ( neighbors->GetId(j) < cellId )
            {
              break;
            }
          }
          if ( j >= numNei )
          {
            edge = VTK_FEATURE_EDGE_VERTEX;
          }
        }
        
        else if ( numNei == 1 && (nei=neighbors->GetId(0)) > cellId ) 
        {
          vtkPolygon::ComputeNormal(inPts,npts,pts,normal);
          Mesh->GetCellPoints(nei,numNeiPts,neiPts);
          vtkPolygon::ComputeNormal(inPts,numNeiPts,neiPts,neiNormal);
          
          if ( this->FeatureEdgeSmoothing &&
               vtkMath::Dot(normal,neiNormal) <= CosFeatureAngle ) 
          {
            edge = VTK_FEATURE_EDGE_VERTEX;
          }
        }
        else // a visited edge; skip rest of analysis
        {
          continue;
        }
        
        if ( edge && Verts[p1].type == VTK_SIMPLE_VERTEX )
        {
          Verts[p1].edges->Reset();
          Verts[p1].edges->InsertNextId(p2);
          Verts[p1].type = edge;
        }
        else if ( (edge && Verts[p1].type == VTK_BOUNDARY_EDGE_VERTEX) ||
                  (edge && Verts[p1].type == VTK_FEATURE_EDGE_VERTEX) ||
                  (!edge && Verts[p1].type == VTK_SIMPLE_VERTEX ) )
        {
          Verts[p1].edges->InsertNextId(p2);
          if ( Verts[p1].type && edge == VTK_BOUNDARY_EDGE_VERTEX )
          {
            Verts[p1].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }
        
        if ( edge && Verts[p2].type == VTK_SIMPLE_VERTEX )
        {
          Verts[p2].edges->Reset();
          Verts[p2].edges->InsertNextId(p1);
          Verts[p2].type = edge;
        }
        else if ( (edge && Verts[p2].type == VTK_BOUNDARY_EDGE_VERTEX ) ||
                  (edge && Verts[p2].type == VTK_FEATURE_EDGE_VERTEX) ||
                  (!edge && Verts[p2].type == VTK_SIMPLE_VERTEX ) )
        {
          Verts[p2].edges->InsertNextId(p1);
          if ( Verts[p2].type && edge == VTK_BOUNDARY_EDGE_VERTEX )
          {
            Verts[p2].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }
      }
    }
    
    inMesh->Delete();
    if (toTris) {toTris->Delete();}
    
    neighbors->Delete();
  }//if strips or polys
  
  //post-process edge vertices to make sure we can smooth them
  for (i=0; i<numPts; i++)
  {
    if ( Verts[i].type == VTK_SIMPLE_VERTEX )
    {
      numSimple++;
    }
    
    else if ( Verts[i].type == VTK_FIXED_VERTEX )
    {
      numFixed++;
    }
    
    else if ( Verts[i].type == VTK_FEATURE_EDGE_VERTEX ||
              Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX )
    { //see how many edges; if two, what the angle is
      
      if ( !this->BoundarySmoothing && 
           Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX )
      {
        Verts[i].type = VTK_FIXED_VERTEX;
        numBEdges++;
      }
      
      else if ( (npts = Verts[i].edges->GetNumberOfIds()) != 2 )
      {
        Verts[i].type = VTK_FIXED_VERTEX;
        numFixed++;
      }
      
      else //check angle between edges
      {
        inPts->GetPoint(Verts[i].edges->GetId(0),x1);
        inPts->GetPoint(i,x2);
        inPts->GetPoint(Verts[i].edges->GetId(1),x3);
        
        for (k=0; k<3; k++)
        {
          l1[k] = x2[k] - x1[k];
          l2[k] = x3[k] - x2[k];
        }
        if ( vtkMath::Normalize(l1) >= 0.0 &&
             vtkMath::Normalize(l2) >= 0.0 &&
             vtkMath::Dot(l1,l2) < CosEdgeAngle)
        {
          numFixed++;
          Verts[i].type = VTK_FIXED_VERTEX;
        }
        else
        {
          if ( Verts[i].type == VTK_FEATURE_EDGE_VERTEX )
          {
            numFEdges++;
          }
          else
          {
            numBEdges++;
          }
        }
      }//if along edge
    }//if edge vertex
  }//for all points
  
  cout<<"Found\n\t" << numSimple << " simple vertices\n\t"
    << numFEdges << " feature edge vertices\n\t"
    << numBEdges << " boundary edge vertices\n\t"
    << numFixed << " fixed vertices\n\t"<<endl;
  cout<<"1:numPts="<<numPts<<endl;
  
  for (i=0; i<numPts; i++) 
  {
    cout<<"Verts["<<i<<"].type="<<VertexType(Verts[i].type)<<endl;
    if(Verts[i].edges != NULL && (npts = Verts[i].edges->GetNumberOfIds()) > 0)
    {
      for (j=0; j<npts; j++)
      {
        cout<<"Verts["<<i<<"].edges->GetId("<<j<<")="<<Verts[i].edges->GetId(j)<<endl;
      }//for all connected points
    }
  }
  
  cout<<"Beginning smoothing iterations..."<<endl;
  cout<<"2:numPts="<<numPts<<endl;
  
  // We've setup the topology...now perform Laplacian smoothing
  //
  newPts = vtkPoints::New();
  newPts->SetNumberOfPoints(numPts);
  
  factor = this->RelaxationFactor;
  for ( maxDist=VTK_DOUBLE_MAX, iterationNumber=0, abortExecute=0; 
        maxDist > conv && iterationNumber < this->NumberOfIterations && !abortExecute;
        iterationNumber++ )
  {
    
    maxDist=0.0;
// 	cout<<"numPts="<<numPts<<endl;
    for (i=0; i<numPts; i++)
    {
//       	cout<<"Verts["<<i<<"].type="<<(int)Verts[i].type<<endl;
      if ( Verts[i].type != VTK_FIXED_VERTEX && Verts[i].edges != NULL &&
           (npts = Verts[i].edges->GetNumberOfIds()) > 0 )
      {
        newPts->GetPoint(i, x); //use current points
        deltaX[0] = deltaX[1] = deltaX[2] = 0.0;
        for (j=0; j<npts; j++)
        {
          newPts->GetPoint(Verts[i].edges->GetId(j), y);
          for (k=0; k<3; k++)
          {
            deltaX[k] += (y[k] - x[k]) / npts;
          }
        }//for all connected points
        
        for (k=0;k<3;k++) 
        {
          xNew[k] = x[k] + factor * deltaX[k];
        }
        
       newPts->SetPoint(i,xNew);
        if ( (dist = vtkMath::Norm(deltaX)) > maxDist )
        {
          maxDist = dist;
        }
      }//if can move point
    }//for all points
  } //for not converged or within iteration count
  
  newPts->Delete();
  
  //free up connectivity storage
  for (i=0; i<numPts; i++)
  {
    if ( Verts[i].edges != NULL ) 
    {
      Verts[i].edges->Delete();
      Verts[i].edges = NULL;
    }
  }
  delete [] Verts;
  
  return 1;
}


CreateSpecialMapping::~CreateSpecialMapping()
{
}


