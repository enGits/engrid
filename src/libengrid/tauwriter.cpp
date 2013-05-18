//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include "tauwriter.h"
#include "guimainwindow.h"
#include <netcdfcpp.h>

//NOTE: on Windows, using "netcdfcpp.h" requires that we retrieve "ncvalues.h"
//directly from ParaView's source code, from the folder "VTK/Utilities/vtknetcdf"


TauWriter::TauWriter()
{
  setFormat("Tau files(*.grid)");
  m_MenuText = "export to Tau format";
}

void TauWriter::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readOutputFileName(file_info.completeBaseName() + ".grid");
    if (isValid()) {
      QString file_name = getFileName();
      NcFile *nc_file = new NcFile(file_name.toAscii(), NcFile::Replace);
      if (!nc_file->is_valid()) {
        EG_ERR_RETURN("unable to open NetCFD file for writing");
      }

      // point coordinates
      size_t Np = m_Grid->GetNumberOfPoints();
      vector<double> x(Np),y(Np),z(Np);
      for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
        vec3_t xv;
        m_Grid->GetPoint(id_node, xv.data());
        x[id_node] = xv[0];
        y[id_node] = xv[1];
        z[id_node] = xv[2];
      }
      NcDim *no_of_points = nc_file->add_dim("no_of_points", Np);
      NcVar *points_xc = nc_file->add_var("points_xc", ncDouble, no_of_points);
      NcVar *points_yc = nc_file->add_var("points_yc", ncDouble, no_of_points);
      NcVar *points_zc = nc_file->add_var("points_zc", ncDouble, no_of_points);
      points_xc->put(&x[0],Np);
      points_yc->put(&y[0],Np);
      points_zc->put(&z[0],Np);

      // boundary faces
      size_t Nbe   = 0;;
      size_t Nquad = 0;
      size_t Ntri  = 0;
      for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
        if (isSurface(id_cell, m_Grid)) {
          ++Nbe;
          if (m_Grid->GetCellType(id_cell) == VTK_TRIANGLE) {
            ++Ntri;
          } else if (m_Grid->GetCellType(id_cell) == VTK_QUAD) {
            ++Nquad;
          } else {
            EG_ERR_RETURN("unsupported boundary element type encountered");
          }
        }
      }
      NcDim *no_of_surfaceelements = nc_file->add_dim("no_of_surfaceelements", Nbe);
      NcVar *boundarymarker_of_surfaces = nc_file->add_var("boundarymarker_of_surfaces", ncInt, no_of_surfaceelements);
      vector<int> bm(Nbe,0), tri(3*Ntri), quad(4*Nquad);
      int i_bm = 0;
      QSet<int> bc_set = getAllBoundaryCodes(m_Grid);
      QVector<int> bcs(bc_set.size());
      qCopy(bc_set.begin(), bc_set.end(), bcs.begin());
      EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
      if (Ntri) {
        NcDim *no_of_surfacetriangles = nc_file->add_dim("no_of_surfacetriangles", Ntri);
        NcDim *points_per_surfacetriangle = nc_file->add_dim("points_per_surfacetriangle", 3);
        NcVar *points_of_surfacetriangles = nc_file->add_var("points_of_surfacetriangles", ncInt, no_of_surfacetriangles, points_per_surfacetriangle);
        int i_tri = 0;
        for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
          if (isSurface(id_cell, m_Grid)) {
            if (m_Grid->GetCellType(id_cell) == VTK_TRIANGLE) {
              for (int i_bc = 0; i_bc < bcs.size(); ++i_bc) {
                if (bcs[i_bc] == cell_code->GetValue(id_cell)) {
                  bm[i_bm] = i_bc+1;
                  break;
                }
              }
              vtkIdType N_pts, *pts;
              m_Grid->GetCellPoints(id_cell, N_pts, pts);
              tri[i_tri + 0] = pts[0];
              tri[i_tri + 1] = pts[1];
              tri[i_tri + 2] = pts[2];
              ++i_bm;
              i_tri += 3;
            }
          }
        }
        points_of_surfacetriangles->put(&tri[0],Ntri,3);
      }
      if (Nquad) {
        NcDim *no_of_surfacequadrilaterals     = nc_file->add_dim("no_of_surfacequadrilaterals",Nquad);
        NcDim *points_per_surfacequadrilateral = nc_file->add_dim("points_per_surfacequadrilateral",4);
        NcVar *points_of_surfacequadrilaterals = nc_file->add_var("points_of_surfacequadrilaterals",ncInt, no_of_surfacequadrilaterals, points_per_surfacequadrilateral);
        int i_quad = 0;
        for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
          if (isSurface(id_cell, m_Grid)) {
            if (m_Grid->GetCellType(id_cell) == VTK_QUAD) {
              for (int i_bc = 0; i_bc < bcs.size(); ++i_bc) {
                if (bcs[i_bc] == cell_code->GetValue(id_cell)) {
                  bm[i_bm] = i_bc+1;
                  break;
                }
              }
              vtkIdType N_pts, *pts;
              m_Grid->GetCellPoints(id_cell, N_pts, pts);
              quad[i_quad + 0] = pts[0];
              quad[i_quad + 1] = pts[1];
              quad[i_quad + 2] = pts[2];
              quad[i_quad + 3] = pts[3];
              ++i_bm;
              i_quad += 4;
            }
          }
        }
        points_of_surfacequadrilaterals->put(&quad[0],Nquad,4);
      }
      boundarymarker_of_surfaces->put(&bm[0],Nbe);
      NcDim *no_of_markers = nc_file->add_dim("no_of_markers", bcs.size());


      int Ntet = 0;
      int Npri = 0;
      int Nhex = 0;
      int Npyr = 0;
      for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
        if (m_Grid->GetCellType(id_cell) == VTK_TETRA)      ++Ntet;
        if (m_Grid->GetCellType(id_cell) == VTK_PYRAMID)    ++Npyr;
        if (m_Grid->GetCellType(id_cell) == VTK_WEDGE)      ++Npri;
        if (m_Grid->GetCellType(id_cell) == VTK_HEXAHEDRON) ++Nhex;
      }

      vector<int> tet(Ntet*4),pyr(Npyr*5),pri(Npri*6),hex(Nhex*8);
      int i_tet = 0;
      int i_pyr = 0;
      int i_pri = 0;
      int i_hex = 0;
      if (Ntet) {
        NcDim *no_of_tetraeders = nc_file->add_dim("no_of_tetraeders",Ntet);
        NcDim *points_per_tetraeder = nc_file->add_dim("points_per_tetraeder",4);
        NcVar *points_of_tetraeders = nc_file->add_var("points_of_tetraeders",ncInt, no_of_tetraeders, points_per_tetraeder);
        for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
          if (m_Grid->GetCellType(id_cell) == VTK_TETRA) {
            vtkIdType *pts, N_pts;
            m_Grid->GetCellPoints(id_cell, N_pts, pts);
            tet[i_tet + 0] = pts[0];
            tet[i_tet + 1] = pts[1];
            tet[i_tet + 2] = pts[2];
            tet[i_tet + 3] = pts[3];
            i_tet += 4;
          }
        }
        points_of_tetraeders->put(&tet[0],Ntet,4);
      }
      if (Npyr) {
        NcDim *no_of_pyramids = nc_file->add_dim("no_of_pyramids",Npyr);
        NcDim *points_per_pyramid = nc_file->add_dim("points_per_pyramid",5);
        NcVar *points_of_pyramids = nc_file->add_var("points_of_pyramids",ncInt, no_of_pyramids, points_per_pyramid);
        for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
          if (m_Grid->GetCellType(id_cell) == VTK_PYRAMID) {
            vtkIdType *pts, N_pts;
            m_Grid->GetCellPoints(id_cell, N_pts, pts);
            pyr[i_pyr + 0] = pts[0];
            pyr[i_pyr + 1] = pts[1];
            pyr[i_pyr + 2] = pts[2];
            pyr[i_pyr + 3] = pts[3];
            pyr[i_pyr + 4] = pts[4];
            i_pyr += 5;
          }
        }
        points_of_pyramids->put(&pyr[0],Npyr,5);
      }
      if (Npri) {
        NcDim *no_of_prisms = nc_file->add_dim("no_of_prisms",Npri);
        NcDim *points_per_prism = nc_file->add_dim("points_per_prism",6);
        NcVar *points_of_prisms = nc_file->add_var("points_of_prisms",ncInt, no_of_prisms,  points_per_prism);
        for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
          if (m_Grid->GetCellType(id_cell) == VTK_WEDGE) {
            vtkIdType *pts, N_pts;
            m_Grid->GetCellPoints(id_cell, N_pts, pts);
            pri[i_pri + 0] = pts[0];
            pri[i_pri + 1] = pts[2];
            pri[i_pri + 2] = pts[1];
            pri[i_pri + 3] = pts[3];
            pri[i_pri + 4] = pts[5];
            pri[i_pri + 5] = pts[4];
            i_pri += 6;
          }
        }
        points_of_prisms->put(&pri[0],Npri,6);
      }
      if (Nhex) {
        NcDim *no_of_hexaeders = nc_file->add_dim("no_of_hexaeders",Nhex);
        NcDim *points_per_hexaeder = nc_file->add_dim("points_per_hexaeder",8);
        NcVar *points_of_hexaeders = nc_file->add_var("points_of_hexaeders",ncInt, no_of_hexaeders, points_per_hexaeder);
        for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
          if (m_Grid->GetCellType(id_cell) == VTK_HEXAHEDRON) {
            vtkIdType *pts, N_pts;
            m_Grid->GetCellPoints(id_cell, N_pts, pts);
            hex[i_hex + 0] = pts[4];
            hex[i_hex + 1] = pts[7];
            hex[i_hex + 2] = pts[6];
            hex[i_hex + 3] = pts[5];
            hex[i_hex + 4] = pts[0];
            hex[i_hex + 5] = pts[3];
            hex[i_hex + 6] = pts[2];
            hex[i_hex + 7] = pts[1];
            i_hex += 8;
          }
        }
        points_of_hexaeders->put(&hex[0],Nhex,8);
      }

      delete nc_file;
    }
  } catch (Error err) {
    err.display();
  }
}

