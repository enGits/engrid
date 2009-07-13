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
#include <cstdlib>

using namespace std;

OpenFOAMTools::OpenFOAMTools(QObject *parent) : QObject( parent )
{
  m_SolverProcess = new QProcess(this);
  m_ToolsProcess = new QProcess(this);

  connect(m_SolverProcess, SIGNAL(error(QProcess::ProcessError)),        this, SLOT(errorHandler(QProcess::ProcessError)));
  connect(m_SolverProcess, SIGNAL(finished(int, QProcess::ExitStatus)),  this, SLOT(finishedHandler(int, QProcess::ExitStatus)));
  connect(m_SolverProcess, SIGNAL(readyReadStandardError()),             this, SLOT(readFromStderr()));
  connect(m_SolverProcess, SIGNAL(readyReadStandardOutput()),            this, SLOT(readFromStdout()));
  connect(m_SolverProcess, SIGNAL(started()),                            this, SLOT(startedHandler()));
  connect(m_SolverProcess, SIGNAL(stateChanged(QProcess::ProcessState)), this, SLOT(stateChangedHandler(QProcess::ProcessState)));

  m_SolverBinary = "";
  m_WorkingDirectory = "";
  m_HostFile = "hostfile.txt";

  QSettings *settings = GuiMainWindow::pointer()->settings();
  if (settings->contains("openfoam_directory")) {
    m_OpenFoamPath = settings->value("openfoam_directory").toString();
  } else {
    m_OpenFoamPath = getenv("HOME");
    m_OpenFoamPath += "/OpenFOAM/OpenFOAM-1.5";
    settings->setValue("openfoam_directory", m_OpenFoamPath);
  }
  if (settings->contains("openfoam_architecture")) {
    m_OpenFoamArch = settings->value("openfoam_architecture").toString();
  } else {
    m_OpenFoamArch = "linux64GccDPOpt";
    settings->setValue("openfoam_architecture", m_OpenFoamArch);
  }

}

OpenFOAMTools::~OpenFOAMTools()
{
  this->stopSolverProcess();
}

void OpenFOAMTools::writeMpiParameters()
{
  QString hostfile_text = GuiMainWindow::pointer()->getXmlSection( "solver/general/host_weight_list" );
  
  QVector <QString> host;
  QVector <QString> weight;
  
  QStringList host_weight_list = hostfile_text.split(",");
  foreach(QString host_weight, host_weight_list) {
    if(!host_weight.isEmpty()){
      QStringList values = host_weight.split(":");
      qWarning()<<"values="<<values;
      host.push_back(values[0].trimmed());
      weight.push_back(values[1].trimmed());
    }
  }
  
  // create the hostfile
  QFileInfo fileinfo( m_WorkingDirectory + "/" + m_HostFile );
  QFile hostfile( fileinfo.filePath() );
  if (!hostfile.open(QIODevice::WriteOnly | QIODevice::Text)) {
    try {
      EG_ERR_RETURN( "ERROR: Failed to open file " + fileinfo.filePath() );
    } catch ( Error err ) {
      err.display();
    }
  }
  QTextStream out( &hostfile );
  for(int i = 0; i < host.size(); i++) {
    out << host[i] << endl;
  }
  hostfile.close();
  
}

int OpenFOAMTools::getArguments()
{
  // resest command-line
  m_Program = "";
  m_Arguments.clear();

  // get binary name
  int solver_type = GuiMainWindow::pointer()->getXmlSection("solver/general/solver_type").toInt();

  QFileInfo solvers_fileinfo;
  solvers_fileinfo.setFile(":/resources/solvers/solvers.txt");
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

  m_SolverBinary = m_OpenFoamPath + "/applications/bin/" + m_OpenFoamArch + "/" + binary;

  if ( m_WorkingDirectory.isEmpty() ) {
    m_WorkingDirectory = QFileDialog::getExistingDirectory( NULL, "select case directory", GuiMainWindow::getCwd() );
    if ( !m_WorkingDirectory.isNull() ) {
      GuiMainWindow::setCwd( QFileInfo( m_WorkingDirectory ).absolutePath() );
    } else {
      return(-1);
    }
  }

  // get number of processes
  m_NumProcesses = GuiMainWindow::pointer()->getXmlSection("solver/general/num_processes");

  // create hostfile + decomposeParDict + get necessary parameters
  writeMpiParameters();
  
  // set working directory of the process
  m_SolverProcess->setWorkingDirectory(m_WorkingDirectory);

  return(0);
}

//=====================================
// Main slots
//=====================================

void OpenFOAMTools::runSolver()
{
  runDecomposePar();
  if (getArguments() < 0) {
    return;
  }
  this->stopSolverProcess();
  if (m_NumProcesses.toInt() <= 1) {
    m_Program = m_SolverBinary;
    m_Arguments << "-case" << m_WorkingDirectory;
  } else {
    m_Program = "mpirun";
    m_Arguments << "--hostfile" << m_HostFile << "-np" << m_NumProcesses << m_SolverBinary << "-case" << m_WorkingDirectory << "-parallel";
  }
  m_SolverProcess->start(m_Program, m_Arguments);
}

void OpenFOAMTools::runTool(QString path, QString name, QStringList args)
{
  m_ToolsProcess->start(m_OpenFoamPath + "/" + path + "/" + m_OpenFoamArch + "/" + name, args);
  do {
    m_ToolsProcess->waitForFinished(500);
    QApplication::processEvents();
  } while (m_ToolsProcess->state() == QProcess::Running);
  cout << m_ToolsProcess->readAllStandardOutput().data();
  flush(cout);
}

void OpenFOAMTools::runDecomposePar()
{
  if (getArguments() < 0) {
    return;
  }
  this->stopSolverProcess();
  m_Program = getBinary("applications/bin", "decomposePar");
  m_Arguments << "-force";
  m_SolverProcess->start(m_Program, m_Arguments);
}

void OpenFOAMTools::runPostProcessingTools()
{
  runTool("applications/bin", "reconstructPar");
  runTool("applications/bin", "foamToVTK");
}

void OpenFOAMTools::stopSolverProcess()
{
  m_SolverProcess->kill();
}

//=====================================
// Handlers
//=====================================

void OpenFOAMTools::finishedHandler(int exitCode, QProcess::ExitStatus exitStatus)
{
  qDebug() << "=== solver-process finished with exit-code = " << exitCode << " ===";
}

void OpenFOAMTools::readFromStderr()
{
  cout << m_SolverProcess->readAllStandardError().data();
  flush(cout);
}

void OpenFOAMTools::readFromStdout()
{
  cout << m_SolverProcess->readAllStandardOutput().data();
  flush(cout);
}

void OpenFOAMTools::startedHandler()
{
  qDebug() << "=== started solver-process with PID = " << m_SolverProcess->pid() << "===";
  QString cmd = m_Program;
  foreach(QString arg, m_Arguments) {
    cmd += " " + arg;
  }
  cout << "[" << qPrintable(m_WorkingDirectory) << "]$ " << qPrintable(cmd) << endl;
}

