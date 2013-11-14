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
#include "volumedefinition.h"
#include "guimainwindow.h"

#include <QFileInfo>
#include <QDir>

FoamWriter::FoamWriter()
{
  EG_TYPENAME;
  setFormat("Foam boundary files(boundary)");
  setExtension("");
  m_CreateCellZones = false;//true;
  m_NoDialog = true;
}

void FoamWriter::writePoints(const PolyMesh &poly)
{
  QString filename = m_Path + "points";
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
  QString filename = m_Path + "faces";
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
  QString filename = m_Path + "owner";
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
  QString filename = m_Path + "neighbour";
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
  QString filename = m_Path + "boundary";
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
    QString bc_type = BC.getType();
    if (GuiMainWindow::pointer()->physicalTypeDefined(bc_type)) {
      PhysicalBoundaryCondition PBC = GuiMainWindow::pointer()->getPhysicalBoundaryCondition(bc_type);
      bc_type = PBC.getFoamType();
    }
    if (hasNeighbour(bc)) {
      bc_type = "mappedWall";
    }

    QString neigh_name = bc_name;

    if (bc_type == "mappedWall") {
      bc_name += "_" + m_CurrentVolume;
    }

    f << "    " << bc_name << "\n";
    f << "    {\n";
    f << "        type                 " << bc_type << ";\n";
    f << "        nFaces               " << nFaces << ";\n";
    f << "        startFace            " << startFace << ";\n";
    if (bc_type == "mappedWall") {
      f << "        startFace            " << startFace << ";\n";
      f << "        sampleMode           nearestPatchFace;\n";
      f << "        sampleRegion         " << getNeighbourName(bc) << ";\n";
      f << "        samplePatch          " << neigh_name + "_" + getNeighbourName(bc) << ";\n";
      f << "        offsetMode           uniform;\n";
      f << "        offset               ( 0 0 0 );\n";
    }
  f << "    }\n";
  };
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
}

void FoamWriter::writeCellZones()
{
  QString filename = m_Path + "cellZones";
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
  f << "    class       regIOobject;\n";
  f << "    location    \"constant/polyMesh\";\n";
  f << "    object      cellZones;\n";
  f << "}\n\n";
  f << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
  f << m_CellZoneLimits.size() - 1 << "\n";
  f << "(\n";
  for (int i = 0; i < m_CellZoneNames.size(); ++i) {
    f << m_CellZoneNames[i] << "\n";
    f << "{\n";
    f << "type cellZone;\n";
    f << "cellLabels List<label>\n";
    f << m_CellZoneLimits[i+1] - m_CellZoneLimits[i] << "\n";
    f << "(\n";
    for (int j = m_CellZoneLimits[i]; j < m_CellZoneLimits[i+1]; ++j) {
      f << j << "\n";
    }
    f << ");\n";
    f << "}\n";
  }
  f << ")\n\n";
  f << "// ************************************************************************* //\n\n\n";
}

PolyMesh* FoamWriter::createSinglePolyMesh()
{
  m_CellZoneLimits.clear();
  m_CellZoneNames.clear();
  QList<VolumeDefinition> vols = mainWindow()->getAllVols();
  PolyMesh* poly = NULL;
  m_CellZoneLimits.append(0);
  for (int i = 0; i < vols.size(); ++i) {
    m_CellZoneNames.append(vols[i].getName());
    EG_VTKSP(vtkUnstructuredGrid, vol_grid);
    MeshPartition volume(vols[i].getName());
    volume.setVolumeOrientation();
    volume.extractToVtkGrid(vol_grid);
    volume.setOriginalOrientation();
    if (i == 0) {
      poly = new PolyMesh(vol_grid, false, 0.0, false);
    } else {
      PolyMesh vol_poly(vol_grid, false, 0.0, false);
      poly->merge(&vol_poly);
    }
    m_CellZoneLimits.append(poly->numPolyCells());
  }
  return poly;
}


void FoamWriter::writeSingleVolume()
{
  bool is_valid = false;
  QString file_name;
  if (!m_NoDialog) {
    readOutputDirectory();
    if (isValid()) {
      is_valid = true;
      file_name = getFileName();
    }
  } else {
    is_valid = true;
    file_name = m_FixedFileName;
  }

  if (is_valid) {
    try {
      QString p1 = file_name;
      QString p2 = p1 + "/constant";
      QDir d1(p1);
      QDir d2(p2);
      if (!d1.exists()) {
        EG_BUG;
      }
      if (!d2.exists()) {
        d1.mkdir("constant");
        d2 = QDir(p2);
      }
      d1 = d2;
      p1 = p2;
      p2 = p1 + "/polyMesh";
      d2 = QDir(p2);
      if (!d2.exists()) {
        d1.mkdir("polyMesh");
      }
      m_Path = file_name + "/constant/polyMesh/";
      if (!QDir(m_Path).exists()) {
        EG_BUG;
      }
      PolyMesh* poly = createSinglePolyMesh();
      writePoints(*poly);
      writeFaces(*poly);
      writeOwner(*poly);
      writeNeighbour(*poly);
      writeBoundary(*poly);
      if (m_CellZoneLimits.size() >= 2) {
        writeCellZones();
      }
      delete poly;
    } catch (Error err) {
      err.display();
    }
  }
}

bool FoamWriter::hasNeighbour(int bc)
{
  if (m_Bc2Vol[bc].size() == 2) {
    return true;
  }
  return false;
}

QString FoamWriter::getNeighbourName(int bc)
{
  foreach (QString name, m_Bc2Vol[bc]) {
    if (name != m_CurrentVolume) {
      return name;
    }
  }
  return "unknown";
}

void FoamWriter::writeMultipleVolumes()
{
  bool is_valid = false;
  QString file_name;
  if (!m_NoDialog) {
    readOutputDirectory();
    if (isValid()) {
      is_valid = true;
      file_name = getFileName();
    }
  } else {
    is_valid = true;
    file_name = m_FixedFileName;
  }

  if (is_valid) {
    try {
      QList<VolumeDefinition> vols = mainWindow()->getAllVols();
      QString p1 = file_name;
      QString p2 = p1 + "/constant";
      QDir d1(p1);
      QDir d2(p2);
      if (!d1.exists()) {
        EG_BUG;
      }
      if (!d2.exists()) {
        d1.mkdir("constant");
        d2 = QDir(p2);
      }

      m_Bc2Vol.clear();
      QSet<int> bcs = mainWindow()->getAllBoundaryCodes();
      foreach (int bc, bcs) {
        foreach (VolumeDefinition vol, vols) {
          if (vol.getSign(bc) != 0) {
            m_Bc2Vol[bc].append(vol.getName());
            if (m_Bc2Vol[bc].size() > 2) {
              EG_ERR_RETURN("Boundary condition with more than two volumes found!");
            }
          }
        }
      }

      foreach (VolumeDefinition vol, vols) {
        m_CurrentVolume = vol.getName();

        QString p3 = p2 + "/" + vol.getName();
        QDir d3(p3);
        if (!d3.exists()) {
          d2.mkdir(vol.getName());
          d3 = QDir(p3);
        }
        QString p4 = p3 + "/polyMesh";
        QDir d4(p4);
        if (!d4.exists()) {
          d3.mkdir("polyMesh");
        }
        m_Path = file_name + "/constant/" + vol.getName() + "/polyMesh/";
        if (!QDir(m_Path).exists()) {
          EG_BUG;
        }
        EG_VTKSP(vtkUnstructuredGrid, vol_grid);
        MeshPartition volume(vol.getName());
        volume.setVolumeOrientation();
        volume.extractToVtkGrid(vol_grid);
        PolyMesh poly(vol_grid, false, 0.0, false);
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
}

void FoamWriter::operate()
{
  if (mainWindow()->getAllVols().size() <= 1 || m_CreateCellZones) {
    writeSingleVolume();
  } else {
    writeMultipleVolumes();
  }
}

void FoamWriter::setFixedFileName(QString file_name)
{
  m_NoDialog = true;
  m_FixedFileName = file_name;
}
