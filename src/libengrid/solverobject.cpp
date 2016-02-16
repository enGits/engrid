// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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

#include "solverobject.h"
#include <iostream>

#include "guimainwindow.h"
#include "filetemplate.h"

SolverObject::SolverObject()
{
  m_CaseDir = "";
  m_BufferedFileName = "";
  m_SolverVersion = "1.5";
}

int SolverObject::deleteBetween(int i, QString str1, QString str2)
{
  int i1 = m_Buffer.indexOf(str1, i);
  if (i1 != -1) {
    int i2 = m_Buffer.indexOf(str2, i1);
    if (i2 != -1) {
      m_Buffer = m_Buffer.remove(i1, i2 - i1 + str2.size());
      return i2-i1;
    }
  }
  return 0;
}

void SolverObject::stripBuffer()
{
  while (deleteBetween(0, "/*", "*/")) {};
  while (deleteBetween(0, "//", "\n")) {};
  int i = m_Buffer.indexOf("FoamFile", 0);
  if (i != -1) {
    deleteBetween(i, "{", "}");
  }
  m_Buffer.replace("FoamFile","");
  m_Buffer.replace("{", " ");
  m_Buffer.replace("}", " ");
  m_Buffer.replace("(", " ");
  m_Buffer.replace(")", " ");
  m_Buffer.replace(";", " ");
  m_Buffer = m_Buffer.simplified();
}

void SolverObject::readFile(QString file_name)
{
  file_name = m_CaseDir + "/" + file_name;
  if(m_BufferedFileName != m_CaseDir + "/" + file_name) {
    m_BufferedFileName = file_name;
    QFile file(file_name);
    if (!file.open(QIODevice::ReadOnly)) {
      EG_ERR_RETURN(QString("error loading file:\n") + file_name);
    }
    QTextStream f(&file);
    m_Buffer = "";
    m_Buffer.reserve(file.size());
    while(!f.atEnd())
    {
      m_Buffer += f.readLine() + "\n";
    }
    stripBuffer();
  }
}

void SolverObject::buildFoamMaps()
{
  int num_foam_nodes;
  readFile("constant/polyMesh/points");
  {
    QTextStream f(getBuffer());
    f >> num_foam_nodes;
  }
  {
    readFile("constant/polyMesh/neighbour");
    QTextStream f(getBuffer());
    int num_neigh;
    f >> num_neigh;
    int neigh = 0;
    m_FirstBoundaryFace = 0;
    while (m_FirstBoundaryFace < num_neigh) {
      f >> neigh;
      if (neigh == -1) {
        break;
      }
      ++m_FirstBoundaryFace;
    }
  }

  m_VolToSurfMap.fill(-1, num_foam_nodes);
  readFile("constant/polyMesh/faces");
  int num_surf_nodes = 0;
  int max_node = 0;
  {
    int num_foam_faces;
    QTextStream f(getBuffer());
    f >> num_foam_faces;
    for (int i = 0; i < num_foam_faces; ++i) {
      int num_nodes;
      f >> num_nodes;
      for (int j = 0; j < num_nodes; ++j) {
        int node;
        f >> node;
        max_node = max(node, max_node);
        if (i >= m_FirstBoundaryFace) {
          if (m_VolToSurfMap[node] == -1) {
            m_VolToSurfMap[node] = num_surf_nodes;
            ++num_surf_nodes;
          }
        }
      }
    }
  }
  m_SurfToVolMap.fill(-1, num_surf_nodes);
  for (int i = 0; i < m_VolToSurfMap.size(); ++i) {
    if (m_VolToSurfMap[i] != -1) {
      if (m_VolToSurfMap[i] > m_SurfToVolMap.size()) {
        EG_BUG;
      }
      m_SurfToVolMap[m_VolToSurfMap[i]] = i;
    }
  }
  for (int i = 0; i < m_SurfToVolMap.size(); ++i) {
    if (m_SurfToVolMap[i] == -1) {
      EG_BUG;
    }
  }
}

void SolverObject::setCaseDir(QString case_dir)
{
  m_CaseDir = case_dir;
  GuiMainWindow::pointer()->setXmlSection("openfoam/CaseDir",m_CaseDir);
}

void SolverObject::writeSolverParameters(QString case_dir)
{
  int idx = GuiMainWindow::pointer()->getXmlSection( "solver/general/solver_type" ).toInt();

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

  QStringList page_list = intext.split( "=" );
  QString page = page_list[idx];
  QString title;
  QString section;
  QString binary;
  QVector <QString> files;
  QStringList variable_list = page.split( ";" );
  foreach( QString variable, variable_list ) {
    QStringList name_value = variable.split( ":" );
    if ( name_value[0].trimmed() == "title" ) title = name_value[1].trimmed();
    if ( name_value[0].trimmed() == "section" ) section = name_value[1].trimmed();
    if ( name_value[0].trimmed() == "binary" ) binary = name_value[1].trimmed();
    if ( name_value[0].trimmed() == "files" ) {
      QStringList file_list = name_value[1].split( "," );
      foreach( QString file, file_list ) {
        files.push_back( file.trimmed() );
      }
    }
    setSolverVersion(title.split(' ').last().trimmed());
  }

  for ( int i = 0; i < files.size(); i++ ) {
    FileTemplate file_template( ":/resources/solvers/" + section + "/" + files[i], section );
    QFileInfo fileinfo_destination( case_dir + "/" + files[i] );
    QDir destination_dir = fileinfo_destination.dir();
    QString destination = case_dir + "/" + files[i];
    if ( !destination_dir.mkpath( destination_dir.absolutePath() ) ) {
      EG_ERR_RETURN( "ERROR: Could not create directory \n" + destination_dir.absolutePath() );
    }
    qDebug() << "Writing to " << destination;
    file_template.exportToOpenFOAM( destination );
  }
}

