//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#ifndef GUIBRLCADIMPORTDIALOGUE_H
#define GUIBRLCADIMPORTDIALOGUE_H

#include <QDialog>

#include "engrid.h"

namespace Ui {
class GuiBrlCadImportDialogue;
}

class GuiBrlCadImportDialogue : public QDialog
{
  Q_OBJECT
  
public:

  explicit GuiBrlCadImportDialogue(QWidget *parent = 0);
  ~GuiBrlCadImportDialogue();

  void    prepare(QString file_name);
  bool    hasSelectedObject();
  QString selectedObject();
  double  scanMemory();

private:

  Ui::GuiBrlCadImportDialogue *ui;

};

#endif // GUIBRLCADIMPORTDIALOGUE_H
