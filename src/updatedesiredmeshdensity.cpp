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
#include "updatedesiredmeshdensity.h"

#include <vtkCharArray.h>

UpdateDesiredMeshDensity::UpdateDesiredMeshDensity()
: SurfaceOperation()
{
  EG_TYPENAME;
}

void UpdateDesiredMeshDensity::operate()
{
  static int nStatic_UpdateDesiredMeshDensity;    // Value of nStatic_UpdateDesiredMeshDensity is retained between each function call
  nStatic_UpdateDesiredMeshDensity++;
  cout << "nStatic_UpdateDesiredMeshDensity is " << nStatic_UpdateDesiredMeshDensity << endl;
  
  //define desired mesh density
  cout<<"=== UpdateDesiredMeshDensity START ==="<<endl;
  cout<<"DebugLevel="<<DebugLevel<<endl;

  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  
  setAllSurfaceCells();
  l2g_t nodes = getPartNodes();
  
  UpdateNodeType();
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  EG_VTKDCN(vtkIntArray, node_specified_density, grid, "node_specified_density");
  
  double diff=Convergence_meshdensity+1;
  if(DebugLevel>3) cout<<"before loop: diff="<<diff<<endl;
  bool first = true;
  int iter = 0;
  do {
    if(DebugLevel>2) cout<<"--->diff="<<diff<<endl;
    first = true;
    foreach (vtkIdType node, nodes) {
      if (DebugLevel > 2) {
        cout << "======>" << endl;
      }
      VertexMeshDensity nodeVMD = getVMD(node);
      int idx = VMDvector.indexOf(nodeVMD);
      node_specified_density->SetValue(node, idx);
      if (DebugLevel > 2) {
        cout << "------>idx=" << idx << endl;
      }
      if ( idx != -1) { //specified
        node_meshdensity_desired->SetValue(node, VMDvector[idx].density);
      } else { //unspecified
        double D = DesiredMeshDensity(node);
        if(first) {
          diff=abs(D-node_meshdensity_desired->GetValue(node));
          first=false;
        } else {
          diff=max(abs(D-node_meshdensity_desired->GetValue(node)),diff);
        }
        node_meshdensity_desired->SetValue(node, D);
      }
      if (DebugLevel > 2) {
        cout << "======>" << endl;
      }
    }
    iter++;
  } while (diff > Convergence_meshdensity && !first && iter < MaxiterDensity);// if first=true, it means no new mesh density has been defined (all densities specified)
  cout << "iter=" << iter << endl;
  if (iter >= MaxiterDensity) {
    cout << "WARNING: Desired convergence factor has not been reached!" << endl;
  }
  
  cout<<"=== UpdateDesiredMeshDensity END ==="<<endl;
}
