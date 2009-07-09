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
#include "filetemplate.h"

#include "guimainwindow.h"

#include <QtDebug>
#include <QValidator>

#include <iostream>
using namespace std;

//=======================================
// Only works if a file is already loaded in engrid (.egc necessary)
int fileTemplateTest( int argc, char ** argv )
{
  QApplication app( argc, argv );

  QVector <QString> files;
  files.push_back( ":/resources/openfoam/simpleFoam/system/fvSchemes.template" );
  files.push_back( ":/resources/openfoam/simpleFoam/system/fvSchemes2.template" );

  TemplateDialog super_gui( files, "openfoam/simplefoam/standard/" );
  super_gui.show();

  return app.exec();
}
//=======================================
// Only works if a file is already loaded in engrid (.egc necessary)
int fileTemplateTest()
{
  QVector <QString> files;
  files.push_back( ":/resources/openfoam/simpleFoam/system/fvSchemes.template" );

  TemplateDialog super_gui( files, "openfoam/simplefoam/standard/" );

  return super_gui.exec();
}
//=======================================
QString TemplateLine::getDefaultValue()
{
  if ( m_DefaultValueEgc == "" ) return m_DefaultValueOpenFOAM;
  else return m_DefaultValueEgc;
}

void TemplateLine::print()
{
  qDebug() << "type=" << this->m_Type;
  qDebug() << "name=" << this->m_Name;
  qDebug() << "options=" << this->m_Options;
  qDebug() << "m_DefaultValueEgc=" << this->m_DefaultValueEgc;
  qDebug() << "m_DefaultValueOpenFOAM=" << this->m_DefaultValueOpenFOAM;
  qDebug() << "position=" << this->m_Position;
}
//=======================================

FileTemplate::FileTemplate()
{
}

FileTemplate::FileTemplate( QString filename, QString section )
{
  this->open( filename, section );
}

void FileTemplate::print()
{
  for ( int i = 0; i < m_Lines.size(); i++ ) {
    m_Lines[i].print();
    qDebug();
  }
}

int FileTemplate::open( QString filename, QString section )
{
  m_Section = section;
  qDebug() << "Opening " << filename;
  m_FileInfo.setFile( filename );
  QFile file( m_FileInfo.filePath() );
  if ( !file.exists() ) {
    qDebug() << "ERROR: " << m_FileInfo.filePath() << " not found.";
    EG_BUG;
    return( -1 );
  }
  if ( !file.open( QIODevice::ReadOnly | QIODevice::Text ) ) {
    qDebug() << "ERROR: Failed to open file.";
    EG_BUG;
    return( -1 );
  }
  QTextStream text_stream( &file );
  m_InText = text_stream.readAll();
  file.close();
  process();
  this->print();
  return( 0 );
}

int FileTemplate::saveEgc()
{
  qDebug() << "Saving EGC ... ";
  QString section = m_Section + m_FileInfo.completeBaseName();
  QString contents = this->getContents();
  GuiMainWindow::pointer()->setXmlSection( section, contents );
  return( 0 );
}

int FileTemplate::exportToOpenFOAM( QString filename )
{
  qDebug() << "Saving openFOAM case as " << filename;

  // set contents
  QString section = m_Section + m_FileInfo.completeBaseName();
  QString openfoam_string = GuiMainWindow::pointer()->getXmlSection( section );
  this->setContents( openfoam_string );

  // save
  m_FileInfo.setFile( filename );
  QFile file( m_FileInfo.filePath() );
  if ( !file.open( QIODevice::WriteOnly | QIODevice::Text ) ) {
    qDebug() << "ERROR: Failed to open file.";
    return( -1 );
  }
  QTextStream out( &file );
  m_OutText = m_InText;
  QRegExp regexp( "<<<.*>>>" );
  regexp.setMinimal( true );
  for ( int i = 0; i < m_Lines.size(); i++ ) {
    int idx1 = m_OutText.indexOf( "<<<" );
    int idx2 = m_OutText.indexOf( ">>>" );
    m_OutText.replace( idx1, idx2 - idx1 + 3, m_Lines[i].getDefaultValue() );
  }
  out << m_OutText;
  file.close();
  return( 0 );
}

int FileTemplate::process()
{
  qDebug() << "Processing...";
  m_Lines.clear();
  QStringList L_open = m_InText.split( "<<<" );
  for ( int i = 1; i < L_open.size(); i++ ) {
    QStringList L_close = L_open[i].split( ">>>" );
    QStringList L_elements = L_close[0].split( ":" );
    TemplateLine template_line;
    template_line.m_Type = L_elements[0];
    template_line.m_Name = L_elements[1];
    template_line.m_Options = L_elements[2];
    template_line.m_DefaultValueOpenFOAM = L_elements[3];
    template_line.m_Position = i;
    m_Lines.push_back( template_line );
  }
  return( 0 );
}

QVector <TemplateLine> FileTemplate::getLines()
{
  return( m_Lines );
}

void FileTemplate::setLines( QVector <TemplateLine> lines )
{
  m_Lines = lines;
}

QString FileTemplate::getContents()
{
  QString ret;
  ret += "\n";
  for ( int i = 0; i < m_Lines.size(); i++ ) {
    ret += m_Lines[i].m_Name + " = " + m_Lines[i].getDefaultValue() + ";\n";
  }
  return ret;
}

void FileTemplate::setContents( QString contents )
{
  QStringList L = contents.split( ";" );
  for ( int i = 0; i < L.size() - 1; i++ ) {
    QStringList L_pair = L[i].split( "=" );
    m_Lines[i].m_DefaultValueEgc = L_pair[1].trimmed();
  }
}
//=======================================

TemplateDialog::TemplateDialog( QVector <QString> files, QString section, QWidget *parent ) : QDialog( parent )
{

  this->setWindowTitle( "Template Viewer" );

  QVBoxLayout* mainLayout = new QVBoxLayout( this );
  this->setLayout( mainLayout );

  QPushButton* openButton = new QPushButton( "Open...", this );
  QPushButton* saveButton = new QPushButton( "Save", this );
  QPushButton* saveAsButton = new QPushButton( "Save as...", this );
  connect( openButton, SIGNAL( clicked() ), this, SLOT( open() ) );
  connect( saveButton, SIGNAL( clicked() ), this, SLOT( saveEgc() ) );
  connect( saveAsButton, SIGNAL( clicked() ), this, SLOT( saveAs() ) );

  for ( int i = 0; i < files.size(); i++ ) {
    m_FileInfo.push_back( QFileInfo( files[i] ) );
    TemplateFormLayout* box = new TemplateFormLayout( files[i], section );
    mainLayout->addLayout( box );
    m_TemplateFormLayoutVector.push_back( box );
  }

  QHBoxLayout *bottomLayout = new QHBoxLayout;
  bottomLayout->addStretch();
  bottomLayout->addWidget( openButton );
  bottomLayout->addWidget( saveButton );
  bottomLayout->addWidget( saveAsButton );
  mainLayout->addLayout( bottomLayout );
}

void TemplateDialog::saveEgc()
{
  for ( int i = 0; i < m_TemplateFormLayoutVector.size(); i++ ) {
    m_TemplateFormLayoutVector[i]->saveEgc();
  }
}

//=======================================

TemplateFormLayout::TemplateFormLayout( QString filename, QString section, char *name, QWidget *parent ) : QFormLayout( parent )
{
  GuiMainWindow::pointer();
  QFormLayout::setObjectName( name );
  m_FileTemplate.open( filename, section );

  QFileInfo file_info( filename );
  QString openfoam_string = GuiMainWindow::pointer()->getXmlSection( section + file_info.completeBaseName() );
  m_FileTemplate.setContents( openfoam_string );

  m_Lines = m_FileTemplate.getLines();
  for ( int i = 0; i < m_Lines.size(); i++ ) {
    if ( m_Lines[i].m_Type == "ComboBox" ) addComboBox( m_Lines[i] );
    else if ( m_Lines[i].m_Type == "IntLineEdit" ) addIntLineEdit( m_Lines[i] );
    else if ( m_Lines[i].m_Type == "DoubleLineEdit" ) addDoubleLineEdit( m_Lines[i] );
    else if ( m_Lines[i].m_Type == "TextLineEdit" ) addTextLineEdit( m_Lines[i] );
    else if ( m_Lines[i].m_Type == "CheckBox" ) addCheckBox( m_Lines[i] );
    else if ( m_Lines[i].m_Type == "SpinBox" ) addSpinBox( m_Lines[i] );
    else if ( m_Lines[i].m_Type == "DoubleSpinBox" ) addDoubleSpinBox( m_Lines[i] );
    else qDebug() << "Unknown type";
  }
}

void TemplateFormLayout::addComboBox( TemplateLine line )
{
  qDebug() << "Adding a ComboBox...";
  QComboBox* combobox = new QComboBox;
  QStringList description;
  QStringList value;
  QStringList L_open = line.m_Options.split( "(" );
  for ( int i = 1; i < L_open.size(); i++ ) {
    QStringList L_close = L_open[i].split( ")" );
    QStringList L_elements = L_close[0].split( "," );
    description << L_elements[0];
    value << L_elements[1].trimmed();
  }
  int current = value.indexOf( line.getDefaultValue().trimmed() );
  qWarning() << "value=" << value;
  qWarning() << "line.getDefaultValue().trimmed()=" << line.getDefaultValue().trimmed();
  qWarning() << "current=" << current;
  combobox->addItems( description );
  combobox-> setCurrentIndex( current );
  this->addRow( line.m_Name, combobox );
  m_ComboBoxVector.push_back( combobox );
  m_ComboboxValues.push_back( value );
}

void TemplateFormLayout::addIntLineEdit( TemplateLine line )
{
  qDebug() << "Adding a IntLineEdit...";
  QValidator *validator = new QIntValidator( this );
  QLineEdit* int_lineedit = new QLineEdit;
  int_lineedit->setValidator( validator );
  int_lineedit->setText( line.getDefaultValue().trimmed() );
  this->addRow( line.m_Name, int_lineedit );
  m_IntLineEditVector.push_back( int_lineedit );
}

void TemplateFormLayout::addDoubleLineEdit( TemplateLine line )
{
  qDebug() << "Adding a DoubleLineEdit...";
  QValidator *validator = new QDoubleValidator( this );
  QLineEdit* double_lineedit = new QLineEdit;
  double_lineedit->setValidator( validator );
  double_lineedit->setText( line.getDefaultValue().trimmed() );
  this->addRow( line.m_Name, double_lineedit );
  m_DoubleLineEditVector.push_back( double_lineedit );
}

void TemplateFormLayout::addTextLineEdit( TemplateLine line )
{
  qDebug() << "Adding a TextLineEdit...";
  QLineEdit* text_lineedit = new QLineEdit;
  text_lineedit->setText( line.getDefaultValue().trimmed() );
  this->addRow( line.m_Name, text_lineedit );
  m_TextLineEditVector.push_back( text_lineedit );
}

void TemplateFormLayout::addCheckBox( TemplateLine line )
{
  qDebug() << "Adding a CheckBox...";
  QCheckBox* check_box = new QCheckBox;
  QStringList L = line.m_Options.split( "," );
  L[0] = L[0].trimmed();
  L[1] = L[1].trimmed();
  int index = L.indexOf( line.getDefaultValue().trimmed() );
  if ( index == 0 ) check_box->setCheckState( Qt::Checked );
  else check_box->setCheckState( Qt::Unchecked );
  QPair < QString, QString > values;
  values.first = L[0];
  values.second = L[1];
  this->addRow( line.m_Name, check_box );
  m_CheckBoxVector.push_back( check_box );
  m_CheckBoxValues.push_back( values );
}

void TemplateFormLayout::addSpinBox( TemplateLine line )
{
  qDebug() << "Adding a SpinBox...";
  QSpinBox* spin_box = new QSpinBox;
  QStringList L = line.m_Options.split( "," );
  int minimum = L[0].trimmed().toInt();
  int maximum = L[1].trimmed().toInt();
  int step = L[2].trimmed().toInt();
  int value = line.getDefaultValue().trimmed().toInt();
  spin_box->setRange( minimum, maximum );
  spin_box->setSingleStep( step );
  spin_box->setValue( value );
  this->addRow( line.m_Name, spin_box );
  m_SpinBoxVector.push_back( spin_box );
}

void TemplateFormLayout::addDoubleSpinBox( TemplateLine line )
{
  qDebug() << "Adding a DoubleSpinBox...";
  QDoubleSpinBox* double_spin_box = new QDoubleSpinBox;
  QStringList L = line.m_Options.split( "," );
  double minimum = L[0].trimmed().toDouble();
  double maximum = L[1].trimmed().toDouble();
  double step = L[2].trimmed().toDouble();
  int decimals = L[3].trimmed().toInt();
  double value = line.getDefaultValue().trimmed().toDouble();
  double_spin_box->setRange( minimum, maximum );
  double_spin_box->setSingleStep( step );
  double_spin_box->setDecimals( decimals );
  double_spin_box->setValue( value );
  this->addRow( line.m_Name, double_spin_box );
  m_DoubleSpinBoxVector.push_back( double_spin_box );
}

void TemplateFormLayout::saveEgc()
{
  getValues();
  m_FileTemplate.saveEgc();
}

void TemplateFormLayout::getValues()
{
  int combobox_idx = 0;
  int intlineedit_idx = 0;
  int doublelineedit_idx = 0;
  int textlineedit_idx = 0;
  int checkbox_idx = 0;
  int spinbox_idx = 0;
  int doublespinbox_idx = 0;
  for ( int i = 0; i < m_Lines.size(); i++ ) {
    if ( m_Lines[i].m_Type == "ComboBox" ) {
      m_Lines[i].m_DefaultValueEgc = readComboBox( combobox_idx );
      combobox_idx++;
    }
    else if ( m_Lines[i].m_Type == "IntLineEdit" ) {
      m_Lines[i].m_DefaultValueEgc =  readIntLineEdit( intlineedit_idx );
      intlineedit_idx++;
    }
    else if ( m_Lines[i].m_Type == "DoubleLineEdit" ) {
      m_Lines[i].m_DefaultValueEgc =  readDoubleLineEdit( doublelineedit_idx );
      doublelineedit_idx++;
    }
    else if ( m_Lines[i].m_Type == "TextLineEdit" ) {
      m_Lines[i].m_DefaultValueEgc =  readTextLineEdit( textlineedit_idx );
      textlineedit_idx++;
    }
    else if ( m_Lines[i].m_Type == "CheckBox" ) {
      m_Lines[i].m_DefaultValueEgc =  readCheckBox( checkbox_idx );
      checkbox_idx++;
    }
    else if ( m_Lines[i].m_Type == "SpinBox" ) {
      m_Lines[i].m_DefaultValueEgc =  readSpinBox( spinbox_idx );
      spinbox_idx++;
    }
    else if ( m_Lines[i].m_Type == "DoubleSpinBox" ) {
      m_Lines[i].m_DefaultValueEgc =  readDoubleSpinBox( doublespinbox_idx );
      doublespinbox_idx++;
    }
    else qDebug() << "Unknown type";
  }
  m_FileTemplate.setLines( m_Lines );
}

QString TemplateFormLayout::readComboBox( int idx )
{
  int i = m_ComboBoxVector[idx]->currentIndex();
  return m_ComboboxValues[idx][i];
}

QString TemplateFormLayout::readIntLineEdit( int idx )
{
  return m_IntLineEditVector[idx]->text();
}

QString TemplateFormLayout::readDoubleLineEdit( int idx )
{
  return m_DoubleLineEditVector[idx]->text();
}

QString TemplateFormLayout::readTextLineEdit( int idx )
{
  return m_TextLineEditVector[idx]->text();
}

QString TemplateFormLayout::readCheckBox( int idx )
{
  if ( m_CheckBoxVector[idx]->checkState() == Qt::Checked ) return m_CheckBoxValues[idx].first;
  else return m_CheckBoxValues[idx].second;
}

QString TemplateFormLayout::readSpinBox( int idx )
{
  QString ret;
  ret.setNum( m_SpinBoxVector[idx]->value() );
  return ret;
}

QString TemplateFormLayout::readDoubleSpinBox( int idx )
{
  QString ret;
  ret.setNum( m_DoubleSpinBoxVector[idx]->value() );
  return ret;
}
//=======================================
