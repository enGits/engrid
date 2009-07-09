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
#include "multipagewidgetpage.h"

#include <QtDebug>

MultiPageWidgetPage::MultiPageWidgetPage( QVector <QString> files, QString section, QWidget *parent )
    : QWidget( parent )
{
  verticalLayout_scrollArea_Solver = new QVBoxLayout( this );
  verticalLayout_scrollArea_Solver->setObjectName( QString::fromUtf8( "verticalLayout_4" ) );
  verticalLayout_scrollArea_Solver->setContentsMargins( 0, 0, 0, 0 );

  scrollArea_Solver = new QScrollArea( this );
  scrollArea_Solver->setObjectName( QString::fromUtf8( "scrollArea_Solver" ) );
  scrollArea_Solver->setWidgetResizable( true );
  scrollArea_Solver->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
  verticalLayout_scrollArea_Solver->addWidget( scrollArea_Solver );

  scrollAreaWidgetContents = new QWidget();
  scrollAreaWidgetContents->setObjectName( QString::fromUtf8( "scrollAreaWidgetContents" ) );
  scrollAreaWidgetContents->setGeometry( QRect( 0, 0, 499, 129 ) );
  scrollArea_Solver->setWidget( scrollAreaWidgetContents );

//   widget = new QWidget( scrollAreaWidgetContents );
//   widget->setObjectName( QString::fromUtf8( "widget" ) );

  verticalLayout_scrollAreaWidgetContents_Solver = new QVBoxLayout();//( widget );
  verticalLayout_scrollAreaWidgetContents_Solver->setObjectName( QString::fromUtf8( "verticalLayout_scrollAreaWidgetContents_Solver" ) );
  verticalLayout_scrollAreaWidgetContents_Solver->setContentsMargins( 0, 0, 0, 0 );
  scrollAreaWidgetContents->setLayout( verticalLayout_scrollAreaWidgetContents_Solver );
  
  for(int i = 0; i < files.size(); i++) {
    TemplateFormLayout* template_form_layout = new TemplateFormLayout(files[i], section);
    verticalLayout_scrollAreaWidgetContents_Solver->addLayout( template_form_layout );
    m_template_form_layout_vector.push_back(template_form_layout);
  }
}

MultiPageWidgetPage::~MultiPageWidgetPage()
{
}

void MultiPageWidgetPage::saveEgc()
{
  qDebug()<<"SAVING PRIVATE RYAN!";
  for(int i = 0; i < m_template_form_layout_vector.size(); i++) {
    m_template_form_layout_vector[i]->saveEgc();
  }
}
