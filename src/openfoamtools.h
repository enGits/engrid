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
#ifndef OPENFOAMTOOLS_H
#define OPENFOAMTOOLS_H

#include <QObject>
#include <QProcess>

#include "egvtkobject.h"

class OpenFOAMTools : public QObject, public EgVtkObject
{
  Q_OBJECT;

private: // attributes

  QProcess*   m_SolverProcess;
  QProcess*   m_ToolsProcess;
  QString     m_SolverBinary;
  QString     m_StrippedSolverBinary;
  QString     m_WorkingDirectory;
  int         m_NumProcesses;
  QString     m_HostFile;
  QString     m_Program;
  QStringList m_Arguments;
  QString     m_OpenFoamPath;
  QString     m_OpenFoamArch;
  QString     m_MainHost;

private: // methods

  void    writeMpiParameters();
  int     getArguments();
  void    runTool(QString path, QString name, QStringList args = QStringList());
  QString getBinary(QString path, QString name) { return m_OpenFoamPath + "/" + path + "/" + m_OpenFoamArch + "/" + name; };

public:

  OpenFOAMTools(QObject *parent = 0);
  ~OpenFOAMTools();

public: // methods



public slots:

  void runSolver();
  void runDecomposePar();
  void runPostProcessingTools();
  void runImportFluentCase();

  void stopSolverProcess();

  void finishedHandler(int exitCode, QProcess::ExitStatus exitStatus);
  void readFromStderr();
  void readFromStdout();
  void startedHandler();

};

#endif
