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

#ifndef GUICREATEHEXCORE_H
#define GUICREATEHEXCORE_H

#include "dialogoperation.h"
#include "ui_guicreatehexcore.h"
#include "createhexcore.h"

class GuiCreateHexCore : public DialogOperation<Ui::GuiCreateHexCore, Operation>
{

  Q_OBJECT

protected: // attributes

  vec3_t m_X1;
  vec3_t m_X2;
  vec3_t m_X10;
  vec3_t m_X20;
  vec3_t m_Xi;
  vec3_t m_Xi0;

protected: // methods

  virtual void before();
  virtual void operate();

public:

  GuiCreateHexCore();

public slots:

  void toggleExternalMesh(bool external_mesh);

};

#endif // GUICREATEHEXCORE_H
