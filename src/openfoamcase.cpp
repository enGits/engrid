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
#include "openfoamcase.h"

#include "filetemplate.h"
#include "guimainwindow.h"

#include <QtDebug>

#include <iostream>
using namespace std;

OpenFOAMcase::OpenFOAMcase()
{
}

void OpenFOAMcase::writeSolverParameters()
{
  int idx = GuiMainWindow::pointer()->getXmlSection("solver/general/solver_type").toInt();
  
  QFileInfo solvers_fileinfo;
  solvers_fileinfo.setFile( ":/resources/solvers/solvers.txt" );
  QFile file( solvers_fileinfo.filePath() );
  if ( !file.exists() ) {
    qDebug() << "ERROR: " << solvers_fileinfo.filePath() << " not found.";
    EG_BUG;
  }
  if ( !file.open( QIODevice::ReadOnly | QIODevice::Text ) ) {
    qDebug() << "ERROR:  Failed to open file " << solvers_fileinfo.filePath();
    EG_BUG;
  }
  QTextStream text_stream( &file );
  QString intext = text_stream.readAll();
  file.close();
  
  QStringList page_list = intext.split("=");
  QString page = page_list[idx];
  QString title;
  QString section;
  QVector <QString> files;
  QStringList variable_list = page.split(";");
  foreach(QString variable, variable_list) {
    QStringList name_value = variable.split(":");
    if(name_value[0].trimmed()=="title") title = name_value[1].trimmed();
    if(name_value[0].trimmed()=="section") section = name_value[1].trimmed();
    if(name_value[0].trimmed()=="files") {
      QStringList file_list = name_value[1].split(",");
      foreach(QString file, file_list) {
        files.push_back(file.trimmed());
      }
    }
  }
  qDebug()<<"title="<<title;
  qDebug()<<"section="<<section;
  qDebug()<<"files="<<files;
  
  for ( int i = 0; i < files.size(); i++ ) {
    FileTemplate file_template( ":/resources/solvers/" + section + "/" + files[i], section );
    QFileInfo fileinfo_destination(getFileName() + "/" + files[i]);
    QDir destination_dir = fileinfo_destination.dir();
    QString destination = destination_dir.absolutePath() + "/" + fileinfo_destination.completeBaseName();
    if(!destination_dir.mkpath(destination_dir.absolutePath())) {
      EG_ERR_RETURN("ERROR: Could not create directory \n" + destination_dir.absolutePath());
    }
    qDebug()<<"Writing to "<<destination;
    file_template.exportToOpenFOAM(destination);
  }
}

void OpenFOAMcase::createBoundaryFaces()
{
  QVector<int> owner(grid->GetNumberOfCells());
  {
    readFile("constant/polyMesh/owner");
    QTextStream f(getBuffer());
    int num_faces;
    f >> num_faces;
    for (int i = 0; i < num_faces; ++i) {
      int o;
      f >> o;
      if (i >= getFirstBoundaryFace()) {
        owner[i - getFirstBoundaryFace()] = o;
      }
    }
  }
  EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
  m_Faces.resize(grid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    face_t F;
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    F.node.resize(N_pts);
    for (int i = 0; i < N_pts; ++i) {
      F.node[i] = surfToVol(pts[i]);
    }
    F.owner = owner[id_cell];
    F.bc = bc->GetValue(id_cell);
    F.neighbour = -1;
    m_Faces[id_cell] = F;
  }
  qSort(m_Faces);
}

void OpenFOAMcase::rewriteBoundaryFaces()
{
  setCaseDir(getFileName());
  buildMaps();
  QVector<QVector<int> > faces(getFirstBoundaryFace());
  {
    readFile("constant/polyMesh/faces");
    QTextStream f(getBuffer());
    int num_faces;
    f >> num_faces;
    if (grid->GetNumberOfCells() + getFirstBoundaryFace() != num_faces) {
      EG_ERR_RETURN("Current surface mesh does not match the OpenFOAM case.");
    }
    for (int i = 0; i < getFirstBoundaryFace(); ++i) {
      int num_nodes;
      f >> num_nodes;
      faces[i].resize(num_nodes);
      for (int j = 0; j < num_nodes; ++j) {
        f >> faces[i][j];
      }
    }
  }
  {
    QString file_name = getFileName() + "/constant/polyMesh/faces";
    QFile file(file_name);
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
    f << faces.size() + grid->GetNumberOfCells() << "\n(\n";
    for (int i = 0; i < faces.size(); ++i) {
      f << faces[i].size() << "(";
      for (int j = 0; j < faces[i].size(); ++j) {
        f << faces[i][j];
        if (j == faces[i].size() - 1) {
          f << ")\n";
        } else {
          f << " ";
        }
      }
    }
    createBoundaryFaces();
    foreach (face_t F, m_Faces) {
      f << F.node.size() << "(";
      for (int i = 0; i < F.node.size(); ++i) {
        f << F.node[i];
        if (i == F.node.size()-1) {
          f << ")\n";
        } else {
          f << " ";
        }
      }
    }
    f << ")\n\n";
    f << "// ************************************************************************* //\n\n\n";
  }
  {
    readFile("constant/polyMesh/owner");
    QTextStream f_in(getBuffer());
    QString file_name = getFileName() + "/constant/polyMesh/owner";
    QFile file(file_name);
    file.open(QIODevice::WriteOnly);
    QTextStream f_out(&file);
    int num_faces;
    f_in >> num_faces;
    f_out << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    f_out << "| =========                 |                                                 |\n";
    f_out << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
    f_out << "|  \\    /   O peration     | Version:  1.5                                   |\n";
    f_out << "|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |\n";
    f_out << "|    \\/     M anipulation  |                                                 |\n";
    f_out << "\\*---------------------------------------------------------------------------*/\n\n";
    f_out << "FoamFile\n";
    f_out << "{\n";
    f_out << "    version     2.0;\n";
    f_out << "    format      ascii;\n";
    f_out << "    class       labelList;\n";
    f_out << "    location    \"constant/polyMesh\";\n";
    f_out << "    object      owner;\n";
    f_out << "}\n\n";
    f_out << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
    f_out << num_faces << "\n(\n";
    for (int i = 0; i < getFirstBoundaryFace(); ++i) {
      int owner;
      f_in >> owner;
      f_out << owner << "\n";
    }
    foreach (face_t F, m_Faces) {
      f_out << F.owner << "\n";
    }
    f_out << ")\n\n";
    f_out << "// ************************************************************************* //\n\n\n";
  }
  m_Path = getFileName() + "/constant/polyMesh/";
  writeBoundary(getFirstBoundaryFace());
}

void OpenFOAMcase::operate()
{
  try {
    if (getFileName() == "") {
      readOutputDirectory();
    }
    if (isValid()) {
      writeSolverParameters();
      rewriteBoundaryFaces();
    }
  } catch (Error err) {
    err.display();
  }
}
