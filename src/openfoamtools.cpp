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
  
  connect( m_Process, SIGNAL(error( QProcess::ProcessError )), this, SLOT(errorHandler()));
  connect( m_Process, SIGNAL(finished ( int , QProcess::ExitStatus )), this, SLOT(finishedHandler()));
  connect( m_Process, SIGNAL(readyReadStandardError()), this, SLOT(readFromStderr()));
  connect( m_Process, SIGNAL(readyReadStandardOutput()), this, SLOT(readFromStdout()));
  connect( m_Process, SIGNAL(started()), this, SLOT(startedHandler()));
  connect( m_Process, SIGNAL(stateChanged( QProcess::ProcessState )), this, SLOT(stateChangedHandler( QProcess::ProcessState )));
}

OpenFOAMTools::~OpenFOAMTools()
{
  this->stopProcesses();
}

void OpenFOAMTools::runSolver()
{
  qDebug()<<"RunSolver";
  m_Process->setWorkingDirectory("/data1/home/mtaverne/Geometries/Testing/tube");
  QString program = "simpleFoam";
  QStringList arguments;
//  arguments << "/data1/home/mtaverne/Geometries/Testing/tube"<<"/data1/home/mtaverne/tmp.txt";
  m_Process->start(program, arguments);
}

void OpenFOAMTools::runFoamToVTK()
{
  qDebug()<<"RunFoamToVTK";
}

void OpenFOAMTools::runDecomposePar()
{
  qDebug()<<"RunDecomposePar";
}

void OpenFOAMTools::runReconstructPar()
{
  qDebug()<<"RunReconstructPar";
}

void OpenFOAMTools::stopProcesses()
{
  qDebug()<<"StopProcesses";
  m_Process->kill();
}

void OpenFOAMTools::errorHandler ( QProcess::ProcessError error )
{

}

void OpenFOAMTools::finishedHandler ( int exitCode, QProcess::ExitStatus exitStatus )
{

}

void OpenFOAMTools::readFromStderr()
{
  qDebug()<<m_Process->readAllStandardError();
}

void OpenFOAMTools::readFromStdout()
{
  qDebug()<<m_Process->readAllStandardOutput();
}

void OpenFOAMTools::startedHandler ()
{

}

void OpenFOAMTools::stateChangedHandler ( QProcess::ProcessState newState )
{

}
