// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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

#ifndef GUIEDGELENGTHSOURCEPIPE_H
#define GUIEDGELENGTHSOURCEPIPE_H

#include "edgelengthsource.h"
#include "ui_guiedgelengthsourcepipe.h"


class GuiEdgeLengthSourcePipe : public GuiEdgeLengthSource<Ui::GuiEdgeLengthSourcePipe>
{
protected: // attributes

  vec3_t m_X1;
  vec3_t m_X2;
  double m_R1;
  double m_R2;
  double m_Length1;
  double m_Length2;


public: // methods

  GuiEdgeLengthSourcePipe();

  virtual bool    read(QString txt);
  virtual QString write();
  virtual void    setDlgFields();
  virtual void    readDlgFields();
  virtual double  edgeLength(vec3_t x);

};

#endif // GUIEDGELENGTHSOURCECONE_H
