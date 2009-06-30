#include "filetemplate.h"
#include <QtDebug>
#include <iostream>

using namespace std;

int fileTemplateTest( int argc, char ** argv ) {
  QApplication app(argc, argv);
  
  GuiTemplateViewer gui_template_viewer("/data1/home/mtaverne/engrid/src/resources/openfoam/simpleFoam/system/fvSchemes");
  gui_template_viewer.show();
  
  return app.exec();
}

void TemplateLine::print() {
  qDebug()<<"type="<<this->type;
  qDebug()<<"name="<<this->name;
  qDebug()<<"options="<<this->options;
  qDebug()<<"position="<<this->position;
}

FileTemplate::FileTemplate() {
}

FileTemplate::FileTemplate(QString filename) {
  this->open(filename);
}

void FileTemplate::print() {
  for(int i = 0; i<m_Lines.size(); i++) {
    m_Lines[i].print();
    qDebug();
  }
}

int FileTemplate::open(QString filename) {
  qDebug()<<"Opening "<<filename;
  m_FileInfo.setFile(filename);
  QFile file(m_FileInfo.filePath());
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    qDebug()<<"ERROR: Failed to open file.";
    return(-1);
  }
  QTextStream text_stream(&file);
  m_InText = text_stream.readAll();
  file.close();
  process();
  return(0);
}

int FileTemplate::process() {
  qDebug()<<"Processing...";
  m_Lines.clear();
  QStringList L_open = m_InText.split("<<<");
  for(int i = 1; i < L_open.size(); i++) {
    QStringList L_close = L_open[i].split(">>>");
    qDebug()<<L_close[0];
    QStringList L_elements = L_close[0].split(":");
    TemplateLine template_line;
    template_line.type = L_elements[0];
    template_line.name = L_elements[1];
    template_line.options = L_elements[2];
    template_line.position = i;
    m_Lines.push_back(template_line);
  }
  return(0);
}

GuiTemplateViewer::GuiTemplateViewer(QWidget *parent) : QDialog(parent) {

}

GuiTemplateViewer::GuiTemplateViewer(QString filename, QWidget *parent) : QDialog(parent) {
  FileTemplate file_template(filename);
  file_template.print();
}

void GuiTemplateViewer::open() {
}

void GuiTemplateViewer::save() {
}

void GuiTemplateViewer::saveAs() {
}
