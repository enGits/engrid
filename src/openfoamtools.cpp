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
  m_NumProcesses = 1;
  m_MainHost = "";
  
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
  
  if(host.size()>=1) {
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
    
  // create the decomposeParDict file
    QFileInfo fileinfo_decomposeParDict( m_WorkingDirectory + "/system/decomposeParDict"  );
    QFile decomposeParDict( fileinfo_decomposeParDict.filePath() );
    if (!decomposeParDict.open(QIODevice::WriteOnly | QIODevice::Text)) {
      try {
        EG_ERR_RETURN( "ERROR: Failed to open file " + fileinfo_decomposeParDict.filePath() );
      } catch ( Error err ) {
        err.display();
      }
    }
    QTextStream out_decomposeParDict( &decomposeParDict );
    out_decomposeParDict << "/*--------------------------------*- C++ -*----------------------------------*\\" << endl;
    out_decomposeParDict << "| =========                 |                                                 |" << endl;
    out_decomposeParDict << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |" << endl;
    out_decomposeParDict << "|  \\    /   O peration     | Version:  1.5                                   |" << endl;
    out_decomposeParDict << "|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |" << endl;
    out_decomposeParDict << "|    \\/     M anipulation  |                                                 |" << endl;
    out_decomposeParDict << "\\*---------------------------------------------------------------------------*/" << endl;
    out_decomposeParDict << "FoamFile" << endl;
    out_decomposeParDict << "{" << endl;
    out_decomposeParDict << "    version     2.0;" << endl;
    out_decomposeParDict << "    format      ascii;" << endl;
    out_decomposeParDict << "    class       dictionary;" << endl;
    out_decomposeParDict << "    object      decomposeParDict;" << endl;
    out_decomposeParDict << "}" << endl;
    out_decomposeParDict << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
    out_decomposeParDict << "" << endl;
    out_decomposeParDict << "numberOfSubdomains "<<weight.size()<<";" << endl;
    out_decomposeParDict << "" << endl;
    out_decomposeParDict << "method          metis;" << endl;
    out_decomposeParDict << "" << endl;
    out_decomposeParDict << "metisCoeffs" << endl;
    out_decomposeParDict << "{" << endl;
    out_decomposeParDict << "    processorWeights " << endl;
    out_decomposeParDict << "    (" << endl;
    for(int i = 0; i < weight.size(); i++) {
      out_decomposeParDict << "        "<<weight[i]<<"" << endl;
    }
    out_decomposeParDict << "    );" << endl;
    out_decomposeParDict << "}" << endl;
    out_decomposeParDict << "" << endl;
    out_decomposeParDict << "distributed     no;" << endl;
    out_decomposeParDict << "" << endl;
    out_decomposeParDict << "roots           " << endl;
    out_decomposeParDict << "(" << endl;
    out_decomposeParDict << ");" << endl;
    out_decomposeParDict << "" << endl;
    out_decomposeParDict << "// ************************************************************************* //" << endl;
    out_decomposeParDict << "" << endl;
    decomposeParDict.close();
    
    m_NumProcesses = weight.size();
    m_MainHost = host[0];
  }
  else {
    m_NumProcesses = 1;
    m_MainHost = "";
  }
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
  if (getArguments() < 0) {
    return;
  }
  this->stopSolverProcess();
  if (m_NumProcesses <= 1) {
    if(m_MainHost.isEmpty()) {
      m_Program = m_SolverBinary;
      m_Arguments << "-case" << m_WorkingDirectory;
    } else {
      m_Program = "ssh";
      m_Arguments << m_MainHost << m_SolverBinary << "-case" << m_WorkingDirectory;
    }
  } else {
    runDecomposePar();
    if(m_SolverProcess->waitForFinished() && m_SolverProcess->exitCode() == 0 ) {
      QString numprocesses_str;
      numprocesses_str.setNum(m_NumProcesses);
      m_Arguments.clear();
      //m_Program = "mpirun";
      //m_Arguments << "--hostfile" << m_HostFile << "-np" << numprocesses_str << m_SolverBinary << "-case" << m_WorkingDirectory << "-parallel";
      m_Program = "ssh";
      m_Arguments << m_MainHost << "mpirun";
      m_Arguments << "--hostfile" << m_WorkingDirectory + "/" + m_HostFile;
      m_Arguments << "-np" << numprocesses_str << m_SolverBinary << "-case" << m_WorkingDirectory << "-parallel";
    }
    else {
      qDebug()<<"ERROR: decomposePar failed.";
      return;
    }
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
  QStringList args;
  args << "-latestTime";
  runTool("applications/bin", "reconstructPar", args);
  runTool("applications/bin", "foamToVTK", args);
}

void OpenFOAMTools::stopSolverProcess()
{
  cout<<"===>void OpenFOAMTools::stopSolverProcess() called"<<endl;
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
