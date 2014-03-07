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

#ifndef EDGELENGTHSOURCE_H
#define EDGELENGTHSOURCE_H

class EdgeLengthSource;
template <class UI> class GuiEdgeLengthSourceDlg;
template <class UI> class GuiEdgeLengthSource;

#include "egvtkobject.h"


class EdgeLengthSource : public EgVtkObject
{

protected: // attributes

  QString m_Name;


public:

  virtual bool    read(QString txt) = 0;
  virtual QString write() = 0;
  virtual void    config() = 0;
  virtual double  edgeLength(vec3_t x) = 0;

  virtual void setName(QString name) = 0;
  QString name() { return m_Name; }

};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <class UI>
class GuiEdgeLengthSourceDlg : public QDialog
{

  friend class GuiEdgeLengthSource<UI>;

protected: // attributes

  UI m_UI; ///< The user interface definition from QtDesigner


public: // methods

  GuiEdgeLengthSourceDlg();

};

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

template <class UI>
GuiEdgeLengthSourceDlg<UI>::GuiEdgeLengthSourceDlg()
{
  m_UI.setupUi(this);
};


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <class UI>
class GuiEdgeLengthSource : public EdgeLengthSource
{

protected: // attributes

  GuiEdgeLengthSourceDlg<UI> *m_Dlg;


protected: // methods

  UI& ui() { return m_Dlg->m_UI; };


public:

  GuiEdgeLengthSource();
  virtual void config();
  virtual void setName(QString name);
  virtual void setDlgFields() = 0;
  virtual void readDlgFields() = 0;

};

// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


template <class UI>
GuiEdgeLengthSource<UI>::GuiEdgeLengthSource()
{
  m_Name = "unnamed";
};

template <class UI>
void GuiEdgeLengthSource<UI>::config()
{
  m_Dlg = new GuiEdgeLengthSourceDlg<UI>;
  setDlgFields();
  if (m_Dlg->exec()) {
    readDlgFields();
  }
  delete m_Dlg;
}

template <class UI>
void GuiEdgeLengthSource<UI>::setName(QString name)
{
  m_Name = name;
}


#endif // EDGELENGTHSOURCE_H
