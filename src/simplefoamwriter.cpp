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
#include "simplefoamwriter.h"

#include <QFileInfo>
#include <QDir>

bool SimpleFoamWriter::face_t::operator<(const face_t &F) const
{
  bool less = false;
  if (bc < F.bc) {
    less = true;
  } else if (bc == F.bc) {
    if (owner < F.owner) {
      less = true;
    } else if (owner == F.owner) {
      if (neighbour < F.neighbour) {
        less = true;
      };
    };
  };
  return less;
};

vec3_t SimpleFoamWriter::face_t::normal(vtkUnstructuredGrid *grid)
{
  if (node.size() < 3) EG_BUG;
  vec3_t xc(0,0,0);
  QVector<vec3_t> x(node.size());
  for (int i = 0; i < node.size(); ++i) {
    grid->GetPoint(node[i],x[i].data());
    xc += x[i];
  };
  xc *= 1.0/node.size();
  vec3_t n(0,0,0);
  for (int i = 0; i < node.size()-1; ++i) {
    vec3_t a = x[i]   - xc; 
    vec3_t b = x[i+1] - xc; 
    n += 0.5*(a.cross(b));
  };
  vec3_t a = x[node.size()-1] - xc; 
  vec3_t b = x[0]             - xc; 
  n += 0.5*(a.cross(b));
  return n;
};


SimpleFoamWriter::SimpleFoamWriter()
{
  setFormat("Foam boundary files(boundary)");
  setExtension("");
}

vtkIdType SimpleFoamWriter::getNeigh(int i_cells, int i_neigh) 
{ 
  int n = c2c[i_cells][i_neigh]; 
  if (n >= 0) return cells[n]; 
  EG_BUG;
  return -1;
};

void SimpleFoamWriter::addFace(face_t F)
{
  if (isVolume(F.neighbour,grid)) {
    if (F.neighbour > F.owner) {
      F.bc = 0;
      lfaces.append(F);
    };
  } else {
    F.bc = bc->GetValue(F.neighbour);
    F.neighbour = -1;
    lfaces.append(F);
  };
};

void SimpleFoamWriter::createFaces()
{
  lfaces.clear();
  EG_VTKDCC(vtkIntArray, cell_code,   grid, "cell_code");
  bc = cell_code;
  eg2of.fill(-1,cells.size());
  int Nvol = 0;
  foreach(int i_cells, cells) {
    vtkIdType *pts;
    vtkIdType  Npts;
    grid->GetCellPoints(cells[i_cells], Npts, pts);
    vtkIdType type_cell = grid->GetCellType(cells[i_cells]);
    vtkIdType id_cell = cells[i_cells];
    
    // tetras
    //
    if (type_cell == VTK_TETRA) {
      eg2of[id_cell] = Nvol++;
      {      
        face_t F(3,id_cell,getNeigh(i_cells,0));
        F.node[0] = pts[2]; F.node[1] = pts[1]; F.node[2] = pts[0];
        addFace(F);
      };
      {
        face_t F(3,id_cell,getNeigh(i_cells,1));
        F.node[0] = pts[0]; F.node[1] = pts[1]; F.node[2] = pts[3];
        addFace(F);
      };
      {
        face_t F(3,id_cell,getNeigh(i_cells,2));
        F.node[0] = pts[0]; F.node[1] = pts[3]; F.node[2] = pts[2];
        addFace(F);
      };
      {
        face_t F(3,id_cell,getNeigh(i_cells,3));
        F.node[0] = pts[1]; F.node[1] = pts[2]; F.node[2] = pts[3];
        addFace(F);
      };
    }
    
    // prisms
    //
    if (type_cell == VTK_WEDGE) {
      eg2of[id_cell] = Nvol++;
      {
        face_t F(3,id_cell,getNeigh(i_cells,0));
        F.node[0] = pts[0]; F.node[1] = pts[1]; F.node[2] = pts[2];
        addFace(F);
      };
      {
        face_t F(3,id_cell,getNeigh(i_cells,1));
        F.node[0] = pts[3]; F.node[1] = pts[5]; F.node[2] = pts[4];
        addFace(F);
      };
      {
        face_t F(4,id_cell,getNeigh(i_cells,2));
        F.node[0] = pts[3]; F.node[1] = pts[4]; F.node[2] = pts[1]; F.node[3] = pts[0];
        addFace(F);
      };
      {
        face_t F(4,id_cell,getNeigh(i_cells,3));
        F.node[0] = pts[1]; F.node[1] = pts[4]; F.node[2] = pts[5]; F.node[3] = pts[2];
        addFace(F);
      };
      {
        face_t F(4,id_cell,getNeigh(i_cells,4));
        F.node[0] = pts[0]; F.node[1] = pts[2]; F.node[2] = pts[5]; F.node[3] = pts[3];
        addFace(F);
      };
    };
    
    // hexes
    //
    if (type_cell == VTK_HEXAHEDRON) {
      eg2of[id_cell] = Nvol++;
      {
        face_t F(4,id_cell,getNeigh(i_cells,0),0);
        F.node[0] = pts[3]; F.node[1] = pts[2]; F.node[2] = pts[1]; F.node[3] = pts[0];
        addFace(F);
      };
      {
        face_t F(4,id_cell,getNeigh(i_cells,1),0);
        F.node[0] = pts[4]; F.node[1] = pts[5]; F.node[2] = pts[6]; F.node[3] = pts[7];
        addFace(F);
      };
      {
        face_t F(4,id_cell,getNeigh(i_cells,2),0);
        F.node[0] = pts[0]; F.node[1] = pts[1]; F.node[2] = pts[5]; F.node[3] = pts[4];
        addFace(F);
      };
      {
        face_t F(4,id_cell,getNeigh(i_cells,3),0);
        F.node[0] = pts[3]; F.node[1] = pts[7]; F.node[2] = pts[6]; F.node[3] = pts[2];
        addFace(F);
      };
      {
        face_t F(4,id_cell,getNeigh(i_cells,4),0);
        F.node[0] = pts[0]; F.node[1] = pts[4]; F.node[2] = pts[7]; F.node[3] = pts[3];
        addFace(F);
      };
      {
        face_t F(4,id_cell,getNeigh(i_cells,5),0);
        F.node[0] = pts[1]; F.node[1] = pts[2]; F.node[2] = pts[6]; F.node[3] = pts[5];
        addFace(F);
      };
    };
  };
  
  faces.resize(lfaces.size());
  qCopy(lfaces.begin(),lfaces.end(),faces.begin());
  qSort(faces);
  
};


void SimpleFoamWriter::writePoints()
{
  QString filename = path + "points";
  QFile file(filename);
  file.open(QIODevice::WriteOnly);
  QTextStream f(&file);
  f << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
  f << "| =========                 |                                                 |\n";
  f << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
  f << "|  \\    /   O peration     | Version:  1.5                                   |\n";
  f << "|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |\n";
  f << "|    \\/     M anipulation  |                                                 |\n";
  f << "\\*---------------------------------------------------------------------------*/\n\n";
  f << "FoamFile\n";
  f << "{\n";
  f << "    version     2.0;\n";
  f << "    format      ascii;\n";
  f << "    class       vectorField;\n";
  f << "    location    \"constant/polyMesh\";\n";
  f << "    object      points;\n";
  f << "}\n\n";
  f << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
  f << nodes.size() << "\n(\n";
  foreach(int i_nodes, nodes) {
    vtkIdType id_node = nodes[i_nodes];
    vec3_t x;
    grid->GetPoint(id_node,x.data());
    f.setRealNumberPrecision(16);
    f << "(" << x[0] << " " << x[1] << " " << x[2] << ")\n";
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
};

void SimpleFoamWriter::writeFaces()
{
  QString filename = path + "faces";
  QFile file(filename);
  file.open(QIODevice::WriteOnly);
  QTextStream f(&file);
  f << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
  f << "| =========                 |                                                 |\n";
  f << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
  f << "|  \\    /   O peration     | Version:  1.5                                   |\n";
  f << "|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |\n";
  f << "|    \\/     M anipulation  |                                                 |\n";
  f << "\\*---------------------------------------------------------------------------*/\n\n";
  f << "FoamFile\n";
  f << "{\n";
  f << "    version     2.0;\n";
  f << "    format      ascii;\n";
  f << "    class       faceList;\n";
  f << "    location    \"constant/polyMesh\";\n";
  f << "    object      faces;\n";
  f << "}\n\n";
  f << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
  f << faces.size() << "\n(\n";
  foreach (face_t F, faces) {
    f << F.node.size() << "(";
    for (int i = 0; i < F.node.size(); ++i) {
      f << F.node[i];
      if (i == F.node.size()-1) {
        f << ")\n";
      } else {
        f << " ";
      };
    };
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
};

void SimpleFoamWriter::writeOwner()
{
  QString filename = path + "owner";
  QFile file(filename);
  file.open(QIODevice::WriteOnly);
  QTextStream f(&file);
  f << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
  f << "| =========                 |                                                 |\n";
  f << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
  f << "|  \\    /   O peration     | Version:  1.5                                   |\n";
  f << "|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |\n";
  f << "|    \\/     M anipulation  |                                                 |\n";
  f << "\\*---------------------------------------------------------------------------*/\n\n";
  f << "FoamFile\n";
  f << "{\n";
  f << "    version     2.0;\n";
  f << "    format      ascii;\n";
  f << "    class       labelList;\n";
  f << "    location    \"constant/polyMesh\";\n";
  f << "    object      owner;\n";
  f << "}\n\n";
  f << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
  f << faces.size() << "\n(\n";
  foreach (face_t F, faces) {
    f << eg2of[F.owner] << "\n";
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
};

void SimpleFoamWriter::writeNeighbour()
{
  QString filename = path + "neighbour";
  QFile file(filename);
  file.open(QIODevice::WriteOnly);
  QTextStream f(&file);
  f << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
  f << "| =========                 |                                                 |\n";
  f << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
  f << "|  \\    /   O peration     | Version:  1.5                                   |\n";
  f << "|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |\n";
  f << "|    \\/     M anipulation  |                                                 |\n";
  f << "\\*---------------------------------------------------------------------------*/\n\n";
  f << "FoamFile\n";
  f << "{\n";
  f << "    version     2.0;\n";
  f << "    format      ascii;\n";
  f << "    class       labelList;\n";
  f << "    location    \"constant/polyMesh\";\n";
  f << "    object      neighbour;\n";
  f << "}\n\n";
  f << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
  f << faces.size() << "\n(\n";
  foreach (face_t F, faces) {
    if (F.neighbour == -1) {
      f << "-1\n";
    } else {
      f << eg2of[F.neighbour] << "\n";
    };
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
};

void SimpleFoamWriter::writeBoundary()
{
  QString filename = path + "boundary";
  QFile file(filename);
  file.open(QIODevice::WriteOnly);
  QTextStream f(&file);
  f << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
  f << "| =========                 |                                                 |\n";
  f << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
  f << "|  \\    /   O peration     | Version:  1.5                                   |\n";
  f << "|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |\n";
  f << "|    \\/     M anipulation  |                                                 |\n";
  f << "\\*---------------------------------------------------------------------------*/\n\n";
  f << "FoamFile\n";
  f << "{\n";
  f << "    version     2.0;\n";
  f << "    format      ascii;\n";
  f << "    class       polyBoundaryMesh;\n";
  f << "    location    \"constant/polyMesh\";\n";
  f << "    object      boundary;\n";
  f << "}\n\n";
  f << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
  QSet<int> bcs;
  foreach (face_t F, faces) {
    if (F.bc != 0) {
      bcs.insert(F.bc);
    };
  };
  f << bcs.size() << "\n(\n";
  QVector<patch_t> patch(bcs.size());
  int N_bc = 0;
  foreach (int bc, bcs) {
    int nFaces = 0;
    int startFace = -1;
    for (int i = 0; i < faces.size(); ++i) {
      if (faces[i].bc == bc) {
        ++nFaces;
        if (startFace == -1) {
          startFace = i;
        };
      };
    };
    if (startFace == -1) {
      EG_BUG;
    }
    patch_t P;
    P.bc = bc;
    P.startFace = startFace;
    P.nFaces = nFaces;
    patch[N_bc] = P;
    ++N_bc;
  }
  qSort(patch);
  foreach (patch_t P, patch) {
    BoundaryCondition BC = getBC(P.bc);
    f << "    " << BC.getName() << "\n";
    f << "    {\n";
    f << "        type        " << BC.getType() << ";\n";
    f << "        nFaces      " << P.nFaces << ";\n";
    f << "        startFace   " << P.startFace << ";\n";
    f << "    }\n";
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
};

void SimpleFoamWriter::operate()
{
  try {
    readOutputDirectory();
    if (isValid()) {
      QString p1 = getFileName();
      QString p2 = p1 + "/constant";
      QDir d1(p1);
      QDir d2(p2);
      if (!d1.exists()) {
        EG_BUG;
      };
      if (!d2.exists()) {
        d1.mkdir("constant");
        d2 = QDir(p2);
      };
      d1 = d2;
      p1 = p2;
      p2 = p1 + "/polyMesh";
      d2 = QDir(p2);
      if (!d2.exists()) {
        d1.mkdir("polyMesh");
      };
      path = getFileName() + "/constant/polyMesh/";
      if (!QDir(path).exists()) {
        EG_BUG;
      };
      createFaces();
      writePoints();
      writeFaces();
      writeOwner();
      writeNeighbour();
      writeBoundary();
    };
  } catch (Error err) {
    err.display();
  };
};
