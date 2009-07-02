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
  
  
  QString str = GuiMainWindow::pointer()->getXmlSection("openfoam/simplefoam/standard");
  qDebug()<<"=======";
  qDebug()<<str;
  qDebug()<<"=======";
  QString str2 = GuiMainWindow::pointer()->getXmlSection("engrid/bc");
  qWarning()<<"=======";
  qWarning()<<str2;
  qWarning()<<"=======";
//   void    setXmlSection(QString name, QString contents);
}
