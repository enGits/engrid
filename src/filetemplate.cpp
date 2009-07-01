#include "filetemplate.h"

#include <QtDebug>
#include <QValidator>

#include <iostream>
using namespace std;

int fileTemplateTest( int argc, char ** argv ) {
  QApplication app(argc, argv);
  
  QVector <QString> files;
  files.push_back("/data1/home/mtaverne/engrid/src/resources/openfoam/simpleFoam/system/fvSchemes.template");
  files.push_back("/data1/home/mtaverne/engrid/src/resources/openfoam/simpleFoam/system/fvSchemes2.template");
  
  SuperGui super_gui(files);
  super_gui.show();
  
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

int FileTemplate::save() {
  this->saveAs(m_FileInfo.completeBaseName());
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

SuperGui::SuperGui(QVector <QString> files, QWidget *parent) : QDialog(parent) {
  
  this->setWindowTitle("Template Viewer");
  
  QVBoxLayout* mainLayout = new QVBoxLayout(this);
  this->setLayout(mainLayout);
  
  QPushButton* openButton = new QPushButton("Open...",this);
  QPushButton* saveButton = new QPushButton("Save",this);
  QPushButton* saveAsButton = new QPushButton("Save as...",this);
  connect(openButton, SIGNAL(clicked()), this, SLOT(open()));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  connect(saveAsButton, SIGNAL(clicked()), this, SLOT(saveAs()));
  
  for(int i = 0; i < files.size(); i++) {
    m_FileInfo.push_back(QFileInfo(files[i]));
    SuperBox* box = new SuperBox(files[i]);
    mainLayout->addLayout(box);
    m_SuperBoxVector.push_back(box);
  }
  
  QHBoxLayout *bottomLayout = new QHBoxLayout;
  bottomLayout->addStretch();
  bottomLayout->addWidget(openButton);
  bottomLayout->addWidget(saveButton);
  bottomLayout->addWidget(saveAsButton);
  mainLayout->addLayout(bottomLayout);
}

SuperBox::SuperBox(QString filename, char *name, QWidget *parent) : QFormLayout(parent) {
  QFormLayout::setObjectName(name);
  m_file_template.open(filename);
  m_file_template.print();
  m_Lines = m_file_template.getLines();
  for(int i = 0; i < m_Lines.size(); i++) {
    if(m_Lines[i].type == "ComboBox") addComboBox(m_Lines[i]);
    else if(m_Lines[i].type == "IntLineEdit") addIntLineEdit(m_Lines[i]);
    else if(m_Lines[i].type == "DoubleLineEdit") addDoubleLineEdit(m_Lines[i]);
    else if(m_Lines[i].type == "TextLineEdit") addTextLineEdit(m_Lines[i]);
    else if(m_Lines[i].type == "CheckBox") addCheckBox(m_Lines[i]);
    else if(m_Lines[i].type == "SpinBox") addSpinBox(m_Lines[i]);
    else if(m_Lines[i].type == "DoubleSpinBox") addDoubleSpinBox(m_Lines[i]);
    else qDebug()<<"Unknown type";
  }
}

void SuperBox::addComboBox(TemplateLine line) {
  qDebug()<<"Adding a ComboBox...";
  QComboBox* combobox = new QComboBox;
  QStringList description;
  QStringList value;
  QStringList L_open = line.options.split("(");
  for(int i = 1; i < L_open.size(); i++) {
    QStringList L_close = L_open[i].split(")");
    QStringList L_elements = L_close[0].split(",");
    description<<L_elements[0];
    value<<L_elements[1];
  }
  combobox->addItems(description);
  this->addRow(line.name, combobox);
  m_ComboBoxVector.push_back(combobox);
  m_ComboboxValues.push_back(value);
}

void SuperBox::addIntLineEdit(TemplateLine line) {
  qDebug()<<"Adding a IntLineEdit...";
  QValidator *validator = new QIntValidator(this);
  QLineEdit* int_lineedit = new QLineEdit;
  int_lineedit->setValidator(validator);
  int_lineedit->setText(line.options.trimmed());
  this->addRow(line.name, int_lineedit);
  m_IntLineEditVector.push_back(int_lineedit);
}

void SuperBox::addDoubleLineEdit(TemplateLine line) {
  qDebug()<<"Adding a DoubleLineEdit...";
  QValidator *validator = new QDoubleValidator(this);
  QLineEdit* double_lineedit = new QLineEdit;
  double_lineedit->setValidator(validator);
  double_lineedit->setText(line.options.trimmed());
  this->addRow(line.name, double_lineedit);
  m_DoubleLineEditVector.push_back(double_lineedit);
}

void SuperBox::addTextLineEdit(TemplateLine line) {
  qDebug()<<"Adding a TextLineEdit...";
  QLineEdit* text_lineedit = new QLineEdit;
  text_lineedit->setText(line.options.trimmed());
  this->addRow(line.name, text_lineedit);
  m_TextLineEditVector.push_back(text_lineedit);
}

void SuperBox::addCheckBox(TemplateLine line) {
  qDebug()<<"Adding a CheckBox...";
  QCheckBox* check_box = new QCheckBox;
  QStringList L = line.options.split(",");
  if(L[0].trimmed() == "checked") check_box->setCheckState(Qt::Checked);
  else check_box->setCheckState(Qt::Unchecked);
  QPair < QString, QString > values;
  values.first = L[1];
  values.second = L[2];
  this->addRow(line.name, check_box);
  m_CheckBoxVector.push_back(check_box);
  m_CheckBoxValues.push_back(values);
}

void SuperBox::addSpinBox(TemplateLine line) {
  qDebug()<<"Adding a SpinBox...";
  QSpinBox* spin_box = new QSpinBox;
  QStringList L = line.options.split(",");
  int minimum = L[0].trimmed().toInt();
  int maximum = L[1].trimmed().toInt();
  int step = L[2].trimmed().toInt();
  int value = L[3].trimmed().toInt();
  spin_box->setRange( minimum, maximum );
  spin_box->setSingleStep( step );
  spin_box->setValue( value );
  this->addRow(line.name, spin_box);
  m_SpinBoxVector.push_back(spin_box);
}

void SuperBox::addDoubleSpinBox(TemplateLine line) {
  qDebug()<<"Adding a DoubleSpinBox...";
  QDoubleSpinBox* double_spin_box = new QDoubleSpinBox;
  QStringList L = line.options.split(",");
  double minimum = L[0].trimmed().toDouble();
  double maximum = L[1].trimmed().toDouble();
  double step = L[2].trimmed().toDouble();
  int decimals = L[3].trimmed().toInt();
  double value = L[4].trimmed().toDouble();
  double_spin_box->setRange( minimum, maximum );
  double_spin_box->setSingleStep( step );
  double_spin_box->setDecimals( decimals );
  double_spin_box->setValue( value );
  this->addRow(line.name, double_spin_box);
  m_DoubleSpinBoxVector.push_back(double_spin_box);
}

void SuperGui::open() {
}

void SuperGui::save() {
  qDebug()<<"Saving...";
  for(int i = 0; i<m_SuperBoxVector.size(); i++) {
    m_SuperBoxVector[i]->save();
  }
}

void SuperGui::saveAs() {
}

void SuperBox::open() {
}

void SuperBox::save() {
  qDebug()<<"Saving...";
  getValues();
  m_file_template.setOutValues(m_OutValues);
  m_file_template.save();
}

void SuperBox::saveAs() {
}

void SuperBox::getValues() {
  int combobox_idx = 0;
  int intlineedit_idx = 0;
  int doublelineedit_idx = 0;
  int textlineedit_idx = 0;
  int checkbox_idx = 0;
  int spinbox_idx = 0;
  int doublespinbox_idx = 0;
  for(int i = 0; i < m_Lines.size(); i++) {
    if(m_Lines[i].type == "ComboBox") {
      m_OutValues<<readComboBox(combobox_idx);
      combobox_idx++;
    }
    else if(m_Lines[i].type == "IntLineEdit") {
      m_OutValues<<readIntLineEdit(intlineedit_idx);
      intlineedit_idx++;
    }
    else if(m_Lines[i].type == "DoubleLineEdit") {
      m_OutValues<<readDoubleLineEdit(doublelineedit_idx);
      doublelineedit_idx++;
    }
    else if(m_Lines[i].type == "TextLineEdit") {
      m_OutValues<<readTextLineEdit(textlineedit_idx);
      textlineedit_idx++;
    }
    else if(m_Lines[i].type == "CheckBox") {
      m_OutValues<<readCheckBox(checkbox_idx);
      checkbox_idx++;
    }
    else if(m_Lines[i].type == "SpinBox") {
      m_OutValues<<readSpinBox(spinbox_idx);
      spinbox_idx++;
    }
    else if(m_Lines[i].type == "DoubleSpinBox") {
      m_OutValues<<readDoubleSpinBox(doublespinbox_idx);
      doublespinbox_idx++;
    }
    else qDebug()<<"Unknown type";
  }
}

QString SuperBox::readComboBox(int idx) {
  int i = m_ComboBoxVector[idx]->currentIndex();
  return m_ComboboxValues[idx][i];
}

QString SuperBox::readIntLineEdit(int idx) {
  return m_IntLineEditVector[idx]->text();
}

QString SuperBox::readDoubleLineEdit(int idx) {
  return m_DoubleLineEditVector[idx]->text();
}

QString SuperBox::readTextLineEdit(int idx) {
  return m_TextLineEditVector[idx]->text();
}

QString SuperBox::readCheckBox(int idx) {
  if(m_CheckBoxVector[idx]->checkState() == Qt::Checked) return m_CheckBoxValues[idx].first;
  else return m_CheckBoxValues[idx].second;
}

QString SuperBox::readSpinBox(int idx) {
  QString ret;
  ret.setNum(m_SpinBoxVector[idx]->value());
  return ret;
}

QString SuperBox::readDoubleSpinBox(int idx) {
  QString ret;
  ret.setNum(m_DoubleSpinBoxVector[idx]->value());
  return ret;
}
