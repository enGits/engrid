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
#include "openfoamtools.h"

#include <QtDebug>

OpenFOAMTools::OpenFOAMTools(QObject *parent)
: QObject(parent)
{
  m_Process = new QProcess(this);
  connect( m_Process, SIGNAL(readyReadStandardOutput()), this, SLOT(readFromStdout()));
}

OpenFOAMTools::~OpenFOAMTools()
{
  this->StopProcesses();
}

void OpenFOAMTools::RunSolver()
{
  qDebug()<<"RunSolver";
  m_Process->setWorkingDirectory("/data1/home/mtaverne/Geometries/Testing/tube");
  QString program = "cd  && ls";
  QStringList arguments;
//  arguments << "/data1/home/mtaverne/Geometries/Testing/tube"<<"/data1/home/mtaverne/tmp.txt";
  m_Process->start(program, arguments);
}

void OpenFOAMTools::RunFoamToVTK()
{
  qDebug()<<"RunFoamToVTK";
}

void OpenFOAMTools::RunDecomposePar()
{
  qDebug()<<"RunDecomposePar";
}

void OpenFOAMTools::RunReconstructPar()
{
  qDebug()<<"RunReconstructPar";
}

void OpenFOAMTools::StopProcesses()
{
  qDebug()<<"StopProcesses";
  m_Process->kill();
}

void OpenFOAMTools::readFromStdout()
{
  qDebug()<<m_Process->readAllStandardOutput();
}
