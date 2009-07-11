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

class OpenFOAMTools : public QObject
{
    Q_OBJECT;

  private:
    QProcess *m_Process;
    QString m_SolverBinary;
    QString m_WorkingDirectory;
    QString m_NumProcesses;
    QString m_HostFile;

  private:
    QString m_Program;
    QStringList m_Arguments;

  public:
    OpenFOAMTools( QObject *parent = 0 );
    ~OpenFOAMTools();

  public:
    int getArguments();

  public slots:
    void runSolver();
    void runFoamToVTK();
    void runDecomposePar();
    void runReconstructPar();
    void stopProcesses();

  public slots:
    void errorHandler( QProcess::ProcessError error );
    void finishedHandler( int exitCode, QProcess::ExitStatus exitStatus );
    void readFromStderr();
    void readFromStdout();
    void startedHandler();
    void stateChangedHandler( QProcess::ProcessState newState );
};

#endif
