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
#ifndef guismoothsurface_H
#define guismoothsurface_H

#include "ui_guismoothsurface.h"
#include "dialogoperation.h"
#include "vertexmeshdensity.h"
#include "settingssheet.h"

#include <vtkPolyDataAlgorithm.h>

class GuiSmoothSurface : public DialogOperation<Ui::GuiSmoothSurface>
{
  
  Q_OBJECT;
  
private slots:
  
  void AddSet();
  void RemoveSet();
  void TestSet();
  void Load();
  void Save();
  void SelectAll_BC();
  void ClearAll_BC();
  void SelectAll_Source();
  void ClearAll_Source();
  
protected: // methods
  
  virtual void before();
  virtual void operate();

private:
  int Nbc;
  SettingsSheet* tableWidget;
public:
  QVector <VertexMeshDensity> GetSet();
  QSettings* local_qset;
  
  /** The currently loaded grid file. */
  QString current_filename;
  
  //  /** The settings file to load. */
  //QString current_settingssheet_name;
  
  int readSettings();
  int writeSettings();
  int DisplayErrorScalars(vtkPolyDataAlgorithm* algo);
  int DisplayErrorVectors(vtkPolyDataAlgorithm* algo);
  
};

#endif
