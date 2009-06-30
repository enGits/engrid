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
  if( !file.exists() ) {
    qDebug()<<"ERROR: "<<m_FileInfo.filePath()<<" not found.";
    return(-1);
  }
  if(!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    qDebug()<<"ERROR: Failed to open file.";
    return(-1);
  }
  QTextStream text_stream(&file);
  m_InText = text_stream.readAll();
  file.close();
  process();
  return(0);
}

int FileTemplate::saveAs(QString filename) {
  qDebug()<<"Saving as "<<filename;
  m_FileInfo.setFile(filename);
  QFile file(m_FileInfo.filePath());
  if(!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    qDebug()<<"ERROR: Failed to open file.";
    return(-1);
  }
  QTextStream out(&file);
  m_OutText = m_InText;
  QRegExp regexp("<<<.*>>>");
  regexp.setMinimal(true);
  for(int i = 0; i < m_OutValues.size(); i++) {
    int idx1 = m_OutText.indexOf("<<<");
    int idx2 = m_OutText.indexOf(">>>");
    m_OutText.replace(idx1,idx2-idx1+3,m_OutValues[i]);
  }
  out<<m_OutText;
  file.close();
  return(0);
}

int FileTemplate::process() {
  qDebug()<<"Processing...";
  m_Lines.clear();
  QStringList L_open = m_InText.split("<<<");
  for(int i = 1; i < L_open.size(); i++) {
    QStringList L_close = L_open[i].split(">>>");
//     qDebug()<<L_close[0];
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

QVector <TemplateLine> FileTemplate::getLines() {
  return(m_Lines);
}

void FileTemplate::setOutValues(QStringList L) {
  m_OutValues = L;
}

GuiTemplateViewer::GuiTemplateViewer(QWidget *parent) : QDialog(parent) {

}

GuiTemplateViewer::GuiTemplateViewer(QString filename, QWidget *parent) : QDialog(parent) {
  
  openButton = new QPushButton("Open...");
  saveButton = new QPushButton("Save");
  saveAsButton = new QPushButton("Save as...");
  
  connect(openButton, SIGNAL(clicked()), this, SLOT(open()));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  connect(saveAsButton, SIGNAL(clicked()), this, SLOT(saveAs()));
  
  QHBoxLayout *bottomLayout = new QHBoxLayout;
  bottomLayout->addStretch();
  bottomLayout->addWidget(openButton);
  bottomLayout->addWidget(saveButton);
  bottomLayout->addWidget(saveAsButton);
  
  formLayout = new QFormLayout;
  mainLayout = new QVBoxLayout;
  
  this->setWindowTitle("Template Viewer");
  
  file_template.open(filename);
  file_template.print();
  m_Lines = file_template.getLines();
  for(int i = 0; i < m_Lines.size(); i++) {
    if(m_Lines[i].type == "ComboBox") addComboBox(m_Lines[i]);
    else qDebug()<<"Unknown type";
  }
  
  mainLayout->addLayout(formLayout);
  mainLayout->addLayout(bottomLayout);
  this->setLayout(mainLayout);
}

void GuiTemplateViewer::addComboBox(TemplateLine line) {
  qDebug()<<"Adding a ComboBox...";
  QComboBox* combobox = new QComboBox;
  QStringList description;
  QStringList value;
  QStringList L_open = line.options.split("(");
  for(int i = 1; i < L_open.size(); i++) {
    QStringList L_close = L_open[i].split(")");
    QStringList L_elements = L_close[0].split(",");
    qDebug()<<L_elements;
    description<<L_elements[0];
    value<<L_elements[1];
  }
  combobox->addItems(description);
  formLayout->addRow(line.name, combobox);
  m_ComboBoxVector.push_back(combobox);
}

void GuiTemplateViewer::getValues() {
  int combobox_idx = 0;
  for(int i = 0; i < m_Lines.size(); i++) {
    if(m_Lines[i].type == "ComboBox") {
      m_OutValues<<readComboBox(combobox_idx);
      combobox_idx++;
    }
    else qDebug()<<"Unknown type";
  }
}

QString GuiTemplateViewer::readComboBox(int idx) {
  return("MOMO");
}

void GuiTemplateViewer::open() {
}

void GuiTemplateViewer::save() {
  qDebug()<<"Saving...";
  getValues();
  file_template.setOutValues(m_OutValues);
  file_template.saveAs("openfoam.txt");
}

void GuiTemplateViewer::saveAs() {
}
