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
  
  for(int i = 0; i<files.size(); i++) {
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
  }
  TemplateDialog super_gui( files );
  super_gui.exec();
  qDebug()<<"GUI DONE";
}
