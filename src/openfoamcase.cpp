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
  files.push_back( "/data1/home/mtaverne/engrid/src/resources/openfoam/simpleFoam/system/fvSchemes.template" );
  files.push_back( "/data1/home/mtaverne/engrid/src/resources/openfoam/simpleFoam/system/fvSchemes2.template" );
  
  FileTemplate file_template(files[0]);
  
  QString section = "openfoam/simplefoam/standard";
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
  }
}
