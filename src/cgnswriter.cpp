#include "cgnswriter.h"

CgnsWriter::CgnsWriter()
{
  setFormat("CGNS files(*.cgns)");
}

void CgnsWriter::writeGrid()
{
  setAllCells();
  eg2cgns.fill(-1, cells.size());

  // create the base node
  if (cg_base_write(fn, "Base", 3, 3, &B)) {
    EG_ERR_RETURN("error creating CGNS base node");
  }

  int Nvcells = 0;
  foreach (vtkIdType id_cell, cells) {
    if (isVolume(id_cell, grid)) {
      ++Nvcells;
    }
  }

  // create the zone node for the main grid
  int size[3];
  size[0] = nodes.size();
  size[1] = Nvcells;
  size[2] = 0;
  if (cg_zone_write(fn,B,"main_grid",size,Unstructured,&Z)) {
    EG_ERR_RETURN("error creating volume zone");
  }

  // create a node for the grid coordinates
  int G;
  if (cg_grid_write(fn,B,Z,"GridCoordinates",&G)) {
    EG_ERR_RETURN("error creating GridCoordinates node");
  }

  // create the arrays for the x,y, and z coordinates
  {
    int C;
    double coord_array[nodes.size()];
    vec3_t x;
    { // x coordinates
      for (int i = 0; i < nodes.size(); ++i) {
        if (nodes[i] != i) {
          EG_ERR_RETURN("unexpected hole in the point list");
        }
        grid->GetPoint(nodes[i], x.data());
        coord_array[i] = x[0];
      }
      if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateX", coord_array, &C)) {
        EG_ERR_RETURN("error writing x-coordinates");
      }
    }
    { // y coordinates
      for (int i = 0; i < nodes.size(); ++i) {
        grid->GetPoint(nodes[i], x.data());
        coord_array[i] = x[1];
      }
      if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateY", coord_array, &C)) {
        EG_ERR_RETURN("error writing y-coordinates");
      }
    }
    { // z coordinates
      for (int i = 0; i < nodes.size(); ++i) {
        grid->GetPoint(nodes[i], x.data());
        coord_array[i] = x[2];
      }
      if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateZ", coord_array, &C)) {
        EG_ERR_RETURN("error writing z-coordinates");
      }
    }
  }

  // count occurrences of different element types
  int ntet = 0;
  int npyr = 0;
  int npri = 0;
  int nhex = 0;
  int ntri = 0;
  int nqua = 0;
  foreach (vtkIdType id_cell, cells) {
    if      (grid->GetCellType(id_cell) == VTK_TETRA)      ++ntet;
    else if (grid->GetCellType(id_cell) == VTK_PYRAMID)    ++npyr;
    else if (grid->GetCellType(id_cell) == VTK_WEDGE)      ++npri;
    else if (grid->GetCellType(id_cell) == VTK_HEXAHEDRON) ++nhex;
    else if (grid->GetCellType(id_cell) == VTK_TRIANGLE)   ++ntri;
    else if (grid->GetCellType(id_cell) == VTK_QUAD)       ++nqua;
  };

  // write element connectivity for tetras
  int start  = 0;
  int end    = 0;
  int i_cgns = 0;
  if (ntet) {
    int S;
    int elements[ntet*4];
    size_t j = 0;
    foreach (vtkIdType id_cell, cells) {
      if (grid->GetCellType(id_cell) == VTK_TETRA) {
        eg2cgns[id_cell] = i_cgns;
        ++i_cgns;
        vtkIdType Npts, *pts;
        grid->GetCellPoints(id_cell, Npts, pts);
        elements[j++] = pts[0]+1;
        elements[j++] = pts[1]+1;
        elements[j++] = pts[2]+1;
        elements[j++] = pts[3]+1;
      }
    }
    start = 1;
    end   = start+ntet-1;
    if (cg_section_write(fn, B, Z, "Tetras", TETRA_4, start, end, 0, elements, &S)) {
      EG_ERR_RETURN("error writing tetras");
    }
  }

  // write element connectivity for pyramids
  if (npyr) {
    int S;
    int elements[npyr*5];
    size_t j = 0;
    foreach (vtkIdType id_cell, cells) {
      if (grid->GetCellType(id_cell) == VTK_PYRAMID) {
        eg2cgns[id_cell] = i_cgns;
        ++i_cgns;
        vtkIdType Npts, *pts;
        grid->GetCellPoints(id_cell, Npts, pts);
        elements[j++] = pts[0]+1;
        elements[j++] = pts[1]+1;
        elements[j++] = pts[2]+1;
        elements[j++] = pts[3]+1;
        elements[j++] = pts[4]+1;
      }
    }
    start = end+1;
    end   = start+npyr-1;
    if (cg_section_write(fn, B, Z, "Pyramids", PYRA_5, start, end, 0, elements, &S)) {
      EG_ERR_RETURN("error writing pyramids");
    }
  }

  // write element connectivity for prisms
  if (npri) {
    int S;
    int elements[npri*6];
    size_t j = 0;
    foreach (vtkIdType id_cell, cells) {
      if (grid->GetCellType(id_cell) == VTK_WEDGE) {
        eg2cgns[id_cell] = i_cgns;
        ++i_cgns;
        vtkIdType Npts, *pts;
        grid->GetCellPoints(id_cell, Npts, pts);
        elements[j++] = pts[3]+1;
        elements[j++] = pts[4]+1;
        elements[j++] = pts[5]+1;
        elements[j++] = pts[0]+1;
        elements[j++] = pts[1]+1;
        elements[j++] = pts[2]+1;
      }
    }
    start = end+1;
    end   = start+npri-1;
    if (cg_section_write(fn, B, Z, "Prisms", PENTA_6, start, end, 0, elements, &S)) {
      EG_ERR_RETURN("error writing prisms");
    }
  }

  // write element connectivity for hexas
  if (nhex) {
    int S;
    int elements[nhex*8];
    size_t j = 0;
    foreach (vtkIdType id_cell, cells) {
      if (grid->GetCellType(id_cell) == VTK_HEXAHEDRON) {
        eg2cgns[id_cell] = i_cgns;
        ++i_cgns;
        vtkIdType Npts, *pts;
        grid->GetCellPoints(id_cell, Npts, pts);
        elements[j++] = pts[0]+1;
        elements[j++] = pts[1]+1;
        elements[j++] = pts[3]+1;
        elements[j++] = pts[2]+1;
        elements[j++] = pts[4]+1;
        elements[j++] = pts[5]+1;
        elements[j++] = pts[7]+1;
        elements[j++] = pts[6]+1;
      }
    }
    start = end+1;
    end   = start+nhex-1;
    if (cg_section_write(fn, B, Z, "Hexes", HEXA_8, start, end, 0, elements, &S)) {
      EG_ERR_RETURN("error writing hexes");
    }
  }

  // write element connectivity for triangles
  if (ntri) {
    int S;
    int elements[ntri*3];
    size_t j = 0;
    foreach (vtkIdType id_cell, cells) {
      if (grid->GetCellType(id_cell) == VTK_TRIANGLE) {
        eg2cgns[id_cell] = i_cgns;
        ++i_cgns;
        vtkIdType Npts, *pts;
        grid->GetCellPoints(id_cell, Npts, pts);
        elements[j++] = pts[2]+1;
        elements[j++] = pts[1]+1;
        elements[j++] = pts[0]+1;
      }
    }
    start = end+1;
    end   = start+ntri-1;
    if (cg_section_write(fn, B, Z, "Triangles", TRI_3, start, end, 0, elements, &S)) {
      EG_ERR_RETURN("error writing triangles");
    }
  }

  // write element connectivity for quads
  if (nqua) {
    int S;
    int elements[nqua*4];
    size_t j = 0;
    foreach (vtkIdType id_cell, cells) {
      if (grid->GetCellType(id_cell) == VTK_QUAD) {
        eg2cgns[id_cell] = i_cgns;
        ++i_cgns;
        vtkIdType Npts, *pts;
        grid->GetCellPoints(id_cell, Npts, pts);
        elements[j++] = pts[3]+1;
        elements[j++] = pts[2]+1;
        elements[j++] = pts[1]+1;
        elements[j++] = pts[0]+1;
      }
    }
    start = end+1;
    end   = start+nqua-1;
    if (cg_section_write(fn, B, Z, "Quads", QUAD_4, start, end, 0, elements, &S)) {
      EG_ERR_RETURN("error writing quads");
    }
  }
}

void CgnsWriter::writeBcs()
{
  EG_VTKDCC(vtkIntArray, cell_code,   grid, "cell_code");
  QSet<int> bcs;
  {
    QVector<vtkIdType> faces;
    getAllSurfaceCells(faces, grid);
    foreach (vtkIdType id_face, faces) {
      bcs.insert(cell_code->GetValue(id_face));
    }
  }
  int BC_cgns;
  foreach (int bc, bcs) {
    QVector<vtkIdType> bc_faces;
    QSet<int> tmp_bcs;
    tmp_bcs.insert(bc);
    getSurfaceCells(tmp_bcs, bc_faces, grid);
    int data[bc_faces.size()];
    for (int i = 0; i < bc_faces.size(); ++i) {
      data[i] = eg2cgns[bc_faces[i]]+1;
    };
    if (cg_boco_write(fn, B, Z, getBC(bc).getName().toAscii().data(), BCTypeNull, ElementList, bc_faces.size(), data, &BC_cgns)) {
      cout << cg_get_error() << endl;
      EG_ERR_RETURN("error writing boundary condition");
    }
  }
}

void CgnsWriter::operate()
{

#ifndef CGNS_SUPPORT
  EG_ERR_RETURN("CGNS support has not been compiled");
#endif

  try {
    readOutputFileName();
    if (isValid()) {
      QString file_name = getFileName();
      if (cg_open(file_name.toAscii().data(), MODE_WRITE, &fn)) {
        EG_ERR_RETURN("error while opening CGNS file  for writing");
      }
      writeGrid();
      writeBcs();
      if (cg_close(fn)) {
        EG_ERR_RETURN("error while closing CGNS file");
      }
    }
  } catch (Error err) {
    err.display();
  }
}

