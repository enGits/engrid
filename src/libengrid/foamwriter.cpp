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
#include "foamwriter.h"

#include <QFileInfo>
#include <QDir>

FoamWriter::FoamWriter()
{
  EG_TYPENAME;
  setFormat("Foam boundary files(boundary)");
  setExtension("");
}

void FoamWriter::writePoints(const PolyMesh &poly)
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
  f << poly.totalNumNodes() << "\n(\n";
  for (int i = 0; i < poly.totalNumNodes(); ++i) {
    vec3_t x = poly.nodeVector(i);
    f.setRealNumberPrecision(16);
    f << "(" << x[0] << " " << x[1] << " " << x[2] << ")\n";
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
}

void FoamWriter::writeFaces(const PolyMesh &poly)
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
  f << poly.numFaces() << "\n(\n";
  for (int i = 0; i < poly.numFaces(); ++i) {
    f << poly.numNodes(i) << "(";
    for (int j = 0; j < poly.numNodes(i); ++j) {
      f << poly.nodeIndex(i,j);
      if (j == poly.numNodes(i) - 1) {
        f << ")\n";
      } else {
        f << " ";
      };
    };
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
}

void FoamWriter::writeOwner(const PolyMesh &poly)
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
  f << poly.numFaces() << "\n(\n";
  for (int i = 0; i < poly.numFaces(); ++i) {
    f << poly.owner(i) << "\n";
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
}

void FoamWriter::writeNeighbour(const PolyMesh &poly)
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
  int N = 0;
  for (int i = 0; i < poly.numFaces(); ++i) {
    if (poly.boundaryCode(i) != 0) break;
    ++N;
  };
  f << N << "\n(\n";
  for (int i = 0; i < N; ++i) {
    f << poly.neighbour(i) << "\n";
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
}

void FoamWriter::writeBoundary(const PolyMesh &poly)
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
  int N = 0;
  for (int i = 0; i < poly.numFaces(); ++i) {
    if (poly.boundaryCode(i) != 0) break;
    ++N;
  };
  f << poly.numBCs() << "\n(\n";
  int i = N;
  while (i < poly.numFaces()) {
    int bc = poly.boundaryCode(i);
    BoundaryCondition BC = getBC(bc);
    int nFaces = 0;
    int startFace = i;
    bool loop = (poly.boundaryCode(i) == bc);
    while (loop) {
      ++nFaces;
      ++i;
      loop = (i < poly.numFaces());
      if (loop) {
        loop = (poly.boundaryCode(i) == bc);
      }
    }
    QString bc_name = BC.getName();
    if (bc_name == "unknown") {
      bc_name.setNum(bc);
      bc_name = "BC_" + bc_name.rightJustified(4, '0');
    }
    f << "    " << bc_name << "\n";
    f << "    {\n";
    f << "        type        " << BC.getType() << ";\n";
    f << "        nFaces      " << nFaces << ";\n";
    f << "        startFace   " << startFace << ";\n";
  f << "    }\n";
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
}

void FoamWriter::operate()
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
      PolyMesh poly(m_Grid);
      writePoints(poly);
      writeFaces(poly);
      writeOwner(poly);
      writeNeighbour(poly);
      writeBoundary(poly);
    }
  } catch (Error err) {
    err.display();
  }
}
