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
#include "openfoamcase.h"

#include "filetemplate.h"
#include "guimainwindow.h"

#include <QtDebug>

#include <iostream>
using namespace std;

OpenFOAMcase::OpenFOAMcase()
    : IOOperation()
{
}

void OpenFOAMcase::operate()
{
  cout << "OpenFOAMcase::operate()" << endl;
  QVector <QString> files;
  files.push_back( ":/resources/openfoam/simpleFoam/system/fvSchemes.template" );
  files.push_back( ":/resources/openfoam/simpleFoam/system/fvSchemes2.template" );
  
/*  for(int i = 0; i<files.size(); i++) {
    QFileInfo file_info(files[i]);
    FileTemplate file_template(files[i]);
    QString section = "openfoam/simplefoam/standard/"+file_info.completeBaseName();
    QString openfoam_string = GuiMainWindow::pointer()->getXmlSection(section);
    qDebug()<<"=======";
    qDebug()<<openfoam_string;
    qDebug()<<"=======";
    if(openfoam_string == "") {
      QString contents = file_template.getContents();
      GuiMainWindow::pointer()->setXmlSection(section, contents);
    }
    else {
      file_template.setContents(openfoam_string);
      qDebug()<<"=== file_template.print(); START ===";
      file_template.print();
      qDebug()<<"=== file_template.print(); END ===";
    }
  }*/
  
  TemplateDialog super_gui( files );
  super_gui.exec();
  qDebug()<<"GUI DONE";
  
  for(int i = 0; i<files.size(); i++) {
    QFileInfo file_info(files[i]);
    FileTemplate file_template(files[i]);
    file_template.save_of();
  }
  
}
