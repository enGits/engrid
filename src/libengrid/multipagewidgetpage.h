// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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
#ifndef MULTIPAGEWIDGETPAGE_H
#define MULTIPAGEWIDGETPAGE_H

#include <QWidget>

#include <QScrollArea>
#include <QVBoxLayout>

#include "filetemplate.h"

class MultiPageWidgetPage : public QWidget
{
    Q_OBJECT
  private:
    QScrollArea *scrollArea_Solver;
    QVBoxLayout *verticalLayout_scrollAreaWidgetContents_Solver;
    QVBoxLayout *verticalLayout_scrollArea_Solver;
    QWidget *scrollAreaWidgetContents;
    QVector <TemplateFormLayout*> m_template_form_layout_vector;
  
  public:
    MultiPageWidgetPage( QVector <QString> files, QString section, QWidget *parent = 0 );
    void saveEgc();
    ~MultiPageWidgetPage();
};

#endif
