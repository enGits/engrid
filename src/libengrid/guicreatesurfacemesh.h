// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include "surfaceoperation.h"
#include "edgelengthsourcemanager.h"

#include <vtkPolyDataAlgorithm.h>

class GuiCreateSurfaceMesh : public DialogOperation<Ui::GuiCreateSurfaceMesh, SurfaceOperation>
{
  
  Q_OBJECT;
  
private slots:
  
  void SelectAll_BC();
  void ClearAll_BC();
  
protected: // methods
  
  virtual void operate();

private:

  int Nbc;
  EdgeLengthSourceManager m_ELSManager;
  int m_NumRows;
  int m_NumCols;
  QVector<QVector<QString> > m_Table;
  void setTextFromTable();
  void getTableFromText();
  
public:

  GuiCreateSurfaceMesh();
  
  QString current_filename;
  
public slots:

  int  readSettings();
  int  writeSettings();
  void read()      { m_ELSManager.read(); }
  void write()     { m_ELSManager.write(); }
  void edit()      { m_ELSManager.edit(); }
  void remove()    { m_ELSManager.remove(); }
  void addSphere() { m_ELSManager.addSphere(); }
  void addCone()   { m_ELSManager.addCone(); }
  void addBox()    { m_ELSManager.addBox(); }
  void addPipe()   { m_ELSManager.addPipe(); }

};

#endif
