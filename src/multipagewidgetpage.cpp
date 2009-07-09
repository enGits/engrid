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

MultiPageWidgetPage::MultiPageWidgetPage( QWidget *parent )
    : QWidget( parent )
{
  verticalLayout_scrollArea_Solver = new QVBoxLayout( this );
  verticalLayout_scrollArea_Solver->setObjectName( QString::fromUtf8( "verticalLayout_4" ) );
  verticalLayout_scrollArea_Solver->setContentsMargins( 0, 0, 0, 0 );

  scrollArea_Solver = new QScrollArea( this );
  scrollArea_Solver->setObjectName( QString::fromUtf8( "scrollArea_Solver" ) );
  scrollArea_Solver->setWidgetResizable( true );
  verticalLayout_scrollArea_Solver->addWidget( scrollArea_Solver );

  scrollAreaWidgetContents = new QWidget();
  scrollAreaWidgetContents->setObjectName( QString::fromUtf8( "scrollAreaWidgetContents" ) );
  scrollAreaWidgetContents->setGeometry( QRect( 0, 0, 499, 129 ) );
  scrollArea_Solver->setWidget( scrollAreaWidgetContents );

  widget = new QWidget( scrollAreaWidgetContents );
  widget->setObjectName( QString::fromUtf8( "widget" ) );

  verticalLayout_scrollAreaWidgetContents_Solver = new QVBoxLayout( widget );
  verticalLayout_scrollAreaWidgetContents_Solver->setObjectName( QString::fromUtf8( "verticalLayout_scrollAreaWidgetContents_Solver" ) );
  verticalLayout_scrollAreaWidgetContents_Solver->setContentsMargins( 0, 0, 0, 0 );
}

MultiPageWidgetPage::~MultiPageWidgetPage()
{
}
