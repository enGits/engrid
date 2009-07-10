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

#include "guimainwindow.h"

#include <QtDebug>
#include <QFileInfo>
#include <QFileDialog>

#include <iostream>
using namespace std;

OpenFOAMTools::OpenFOAMTools( QObject *parent )
    : QObject( parent )
{
  m_Process = new QProcess( this );

  connect( m_Process, SIGNAL( error( QProcess::ProcessError ) ), this, SLOT( errorHandler( QProcess::ProcessError ) ) );
  connect( m_Process, SIGNAL( finished( int , QProcess::ExitStatus ) ), this, SLOT( finishedHandler( int , QProcess::ExitStatus ) ) );
  connect( m_Process, SIGNAL( readyReadStandardError() ), this, SLOT( readFromStderr() ) );
  connect( m_Process, SIGNAL( readyReadStandardOutput() ), this, SLOT( readFromStdout() ) );
  connect( m_Process, SIGNAL( started() ), this, SLOT( startedHandler() ) );
  connect( m_Process, SIGNAL( stateChanged( QProcess::ProcessState ) ), this, SLOT( stateChangedHandler( QProcess::ProcessState ) ) );

  m_SolverBinary = "";
  m_WorkingDirectory = "";
  m_HostFile = "hostfile.txt";
}

OpenFOAMTools::~OpenFOAMTools()
{
  this->stopProcesses();
}

int OpenFOAMTools::getArguments()
{
  // resest command-line
  m_Program = "";
  m_Arguments.clear();

  // get binary name
  int solver_type = GuiMainWindow::pointer()->getXmlSection( "solver/general/solver_type" ).toInt();

  QFileInfo solvers_fileinfo;
  solvers_fileinfo.setFile( ":/resources/solvers/solvers.txt" );
  QFile file( solvers_fileinfo.filePath() );
  if ( !file.exists() ) {
    qDebug() << "ERROR: " << solvers_fileinfo.filePath() << " not found.";
    EG_BUG;
    return( -1 );
  }
  if ( !file.open( QIODevice::ReadOnly | QIODevice::Text ) ) {
    qDebug() << "ERROR:  Failed to open file " << solvers_fileinfo.filePath();
    EG_BUG;
    return( -1 );
  }
  QTextStream text_stream( &file );
  QString intext = text_stream.readAll();
  file.close();

  QStringList page_list = intext.split( "=" );
  QString page = page_list[solver_type];
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
  }

  m_SolverBinary = binary;

  // get case directory if it is undefined
//   m_WorkingDirectory = GuiMainWindow::getCwd();
  if ( m_WorkingDirectory.isEmpty() ) {
    m_WorkingDirectory = QFileDialog::getExistingDirectory( NULL, "write OpenFOAM mesh", GuiMainWindow::getCwd() );
    if ( !m_WorkingDirectory.isNull() ) {
      GuiMainWindow::setCwd( QFileInfo( m_WorkingDirectory ).absolutePath() );
    }
    else {
      return( -1 );
    }
  }

  // get number of processes
  m_NumProcesses = GuiMainWindow::pointer()->getXmlSection( "solver/general/num_processes" );

  // create the hostfile
  QString hostfile_text = GuiMainWindow::pointer()->getXmlSection( "solver/general/hostfile" );

  QFileInfo fileinfo( m_WorkingDirectory + "/" + m_HostFile );
  QFile hostfile( fileinfo.filePath() );
  if ( !hostfile.open( QIODevice::WriteOnly | QIODevice::Text ) ) {
    try {
      EG_ERR_RETURN( "ERROR: Failed to open file " + fileinfo.filePath() );
    }
    catch ( Error err ) {
      err.display();
      return( -1 );
    }
  }
  QTextStream out( &hostfile );
  out << hostfile_text;
  hostfile.close();

  // set working directory of the process
  m_Process->setWorkingDirectory( m_WorkingDirectory );

  return( 0 );
}

//=====================================
// Main slots
//=====================================

void OpenFOAMTools::runSolver()
{
  qDebug() << "=== RunSolver ===";
  if ( getArguments() < 0 ) return;
  this->stopProcesses();

  if ( m_NumProcesses.toInt() <= 1 ) {
    m_Program = m_SolverBinary;
    m_Arguments << "-case" << m_WorkingDirectory;
  }
  else {
    m_Program = "mpirun";
    m_Arguments << "--hostfile" << m_HostFile << "-np" << m_NumProcesses << m_SolverBinary << "-case" << m_WorkingDirectory << "-parallel";
  }
  m_Process->start( m_Program, m_Arguments );
}

void OpenFOAMTools::runFoamToVTK()
{
  qDebug() << "=== RunFoamToVTK ===";
  if ( getArguments() < 0 ) return;
  this->stopProcesses();

  m_Program = "foamToVTK";
  m_Process->start( m_Program, m_Arguments );
}

void OpenFOAMTools::runDecomposePar()
{
  qDebug() << "=== RunDecomposePar ===";
  if ( getArguments() < 0 ) return;
  this->stopProcesses();

  m_Program = "decomposePar";
  m_Process->start( m_Program, m_Arguments );
}

void OpenFOAMTools::runReconstructPar()
{
  qDebug() << "=== RunReconstructPar ===";
  if ( getArguments() < 0 ) return;
  this->stopProcesses();

  m_Program = "reconstructPar";
  m_Process->start( m_Program, m_Arguments );
}

void OpenFOAMTools::stopProcesses()
{
  qDebug() << "=== Stopping processes ===";
  m_Process->kill();
}

//=====================================
// Handlers
//=====================================

void OpenFOAMTools::errorHandler( QProcess::ProcessError error )
{
  qDebug() << "=== A process error occured. ===";
}

void OpenFOAMTools::finishedHandler( int exitCode, QProcess::ExitStatus exitStatus )
{
  qDebug() << "=== Process finished. ===";
}

void OpenFOAMTools::readFromStderr()
{
//   qDebug()<<m_Process->readAllStandardError();
  cout << m_Process->readAllStandardError().data() << endl;
}

void OpenFOAMTools::readFromStdout()
{
//   qDebug()<<m_Process->readAllStandardOutput();
  cout << m_Process->readAllStandardOutput().data() << endl;
}

void OpenFOAMTools::startedHandler()
{
  qDebug() << "=== Started process with PID = " << m_Process->pid() << "===";
  QString cmd = m_Program;
  foreach( QString arg, m_Arguments ) cmd += " " + arg;
  cout << "[" << qPrintable( m_WorkingDirectory ) << "]$ " << qPrintable( cmd ) << endl;
}

void OpenFOAMTools::stateChangedHandler( QProcess::ProcessState newState )
{
  qDebug() << "=== Process state changed. ===";
}
