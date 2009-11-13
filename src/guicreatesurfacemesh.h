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
#ifndef guicreatesurfacemesh_H
#define guicreatesurfacemesh_H

#include "ui_guicreatesurfacemesh.h"
#include "dialogoperation.h"
#include "vertexmeshdensity.h"
#include "settingssheet.h"
#include "surfaceoperation.h"
#include "edgelengthsourcemanager.h"

#include <vtkPolyDataAlgorithm.h>

class GuiCreateSurfaceMesh : public DialogOperation<Ui::GuiCreateSurfaceMesh, SurfaceOperation>
{
  
  Q_OBJECT;
  
private slots:
  
  void AddSet();
  void RemoveSet();
  void TestSet();
  void SelectAll_BC();
  void ClearAll_BC();
  
protected: // methods
  
  virtual void operate();

private:

  int Nbc;
  SettingsSheet* m_tableWidget;
  EdgeLengthSourceManager m_ELSManager;
  
public:

  GuiCreateSurfaceMesh();
  
  QVector <VertexMeshDensity> getSet();
  
  QString current_filename;
  
  int DisplayErrorScalars(vtkPolyDataAlgorithm* algo);
  int DisplayErrorVectors(vtkPolyDataAlgorithm* algo);

public slots:

  int  readSettings();
  int  writeSettings();
  void read()      { m_ELSManager.read(); };
  void write()     {m_ELSManager.write(); };
  void edit()      {m_ELSManager.edit(); };
  void remove()    {m_ELSManager.remove(); };
  void addSphere() {m_ELSManager.addSphere();; };

};

#endif
