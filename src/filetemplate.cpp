#include "filetemplate.h"

#include "guimainwindow.h"

#include <QtDebug>
#include <QValidator>

#include <iostream>
using namespace std;

//=======================================
int fileTemplateTest( int argc, char ** argv )
{
  QApplication app( argc, argv );

  QVector <QString> files;
  files.push_back( "/data1/home/mtaverne/engrid/src/resources/openfoam/simpleFoam/system/fvSchemes.template" );
  files.push_back( "/data1/home/mtaverne/engrid/src/resources/openfoam/simpleFoam/system/fvSchemes2.template" );

  TemplateDialog super_gui( files );
  super_gui.show();

  return app.exec();
}
//=======================================
QString TemplateLine::getDefaultValue()
{
  if(default_value_egc == "") return default_value_of;
  else return default_value_egc;
}

void TemplateLine::print()
{
  qDebug() << "type=" << this->type;
  qDebug() << "name=" << this->name;
  qDebug() << "options=" << this->options;
  qDebug() << "default_value_egc=" << this->default_value_egc;
  qDebug() << "default_value_of=" << this->default_value_of;
  qDebug() << "position=" << this->position;
}
//=======================================

FileTemplate::FileTemplate()
{
}

FileTemplate::FileTemplate( QString filename )
{
  this->open( filename );
}

void FileTemplate::print()
{
  for ( int i = 0; i < m_Lines.size(); i++ ) {
    m_Lines[i].print();
    qDebug();
  }
}

int FileTemplate::open( QString filename )
{
  qDebug() << "Opening " << filename;
  m_FileInfo.setFile( filename );
  QFile file( m_FileInfo.filePath() );
  if ( !file.exists() ) {
    qDebug() << "ERROR: " << m_FileInfo.filePath() << " not found.";
    return( -1 );
  }
  if ( !file.open( QIODevice::ReadOnly | QIODevice::Text ) ) {
    qDebug() << "ERROR: Failed to open file.";
    return( -1 );
  }
  QTextStream text_stream( &file );
  m_InText = text_stream.readAll();
  file.close();
  process();
  return( 0 );
}

int FileTemplate::save_egc()
{
  qDebug() << "Saving EGC ... ";
  QString section = "openfoam/simplefoam/standard/"+m_FileInfo.completeBaseName();
  QString contents = this->getContents();
  GuiMainWindow::pointer()->setXmlSection(section, contents);
  
  return( 0 );
}

int FileTemplate::save_of()
{
  qDebug() << "Saving OF ...";
  
  
  QString section = "openfoam/simplefoam/standard/"+m_FileInfo.completeBaseName();
  QString openfoam_string = GuiMainWindow::pointer()->getXmlSection(section);
  this->setContents(openfoam_string);
  qDebug()<<"=== After reading EGC START ===";
  this->print();
  qDebug()<<"=== After reading EGC END ===";
  
  this->saveAs_of( m_FileInfo.completeBaseName() );
  return( 0 );
}

int FileTemplate::saveAs_of( QString filename )
{
  qDebug() << "Saving as " << filename;
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
    template_line.type = L_elements[0];
    template_line.name = L_elements[1];
    template_line.options = L_elements[2];
    template_line.default_value_of = L_elements[3];
    template_line.position = i;
    m_Lines.push_back( template_line );
  }
  return( 0 );
}

QVector <TemplateLine> FileTemplate::getLines()
{
  return( m_Lines );
}

QString FileTemplate::getContents() {
  QString ret;
  ret += "\n";
  for(int i = 0; i < m_Lines.size(); i++) {
    ret += m_Lines[i].name + " = " + m_Lines[i].getDefaultValue() + ";\n";
  }
  return ret;
}

void FileTemplate::setContents(QString contents) {
  QStringList L = contents.split(";");
  for(int i = 0; i < L.size()-1; i++) {
    QStringList L_pair = L[i].split("=");
    m_Lines[i].default_value_egc = L_pair[1].trimmed();
  }
}
//=======================================

TemplateDialog::TemplateDialog( QVector <QString> files, QWidget *parent ) : QDialog( parent )
{

  this->setWindowTitle( "Template Viewer" );

  QVBoxLayout* mainLayout = new QVBoxLayout( this );
  this->setLayout( mainLayout );

  QPushButton* openButton = new QPushButton( "Open...", this );
  QPushButton* saveButton = new QPushButton( "Save", this );
  QPushButton* saveAsButton = new QPushButton( "Save as...", this );
  connect( openButton, SIGNAL( clicked() ), this, SLOT( open() ) );
  connect( saveButton, SIGNAL( clicked() ), this, SLOT( save() ) );
  connect( saveAsButton, SIGNAL( clicked() ), this, SLOT( saveAs() ) );

  for ( int i = 0; i < files.size(); i++ ) {
    m_FileInfo.push_back( QFileInfo( files[i] ) );
    TemplateFormLayout* box = new TemplateFormLayout( files[i] );
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

void TemplateDialog::open()
{
}

void TemplateDialog::save()
{
  qDebug() << "Saving...";
  for ( int i = 0; i < m_TemplateFormLayoutVector.size(); i++ ) {
    m_TemplateFormLayoutVector[i]->save();
  }
}

void TemplateDialog::saveAs()
{
}

//=======================================

TemplateFormLayout::TemplateFormLayout( QString filename, char *name, QWidget *parent ) : QFormLayout( parent )
{
  GuiMainWindow::pointer();
  QFormLayout::setObjectName( name );
  m_file_template.open( filename );
  qDebug()<<"=== Before reading EGC START ===";
  m_file_template.print();
  qDebug()<<"=== Before reading EGC END ===";
  
  QFileInfo file_info(filename);
  QString section = "openfoam/simplefoam/standard/"+file_info.completeBaseName();
  QString openfoam_string = GuiMainWindow::pointer()->getXmlSection(section);
  m_file_template.setContents(openfoam_string);
  qDebug()<<"=== After reading EGC START ===";
  m_file_template.print();
  qDebug()<<"=== After reading EGC END ===";
    
  m_Lines = m_file_template.getLines();
  for ( int i = 0; i < m_Lines.size(); i++ ) {
    if ( m_Lines[i].type == "ComboBox" ) addComboBox( m_Lines[i] );
    else if ( m_Lines[i].type == "IntLineEdit" ) addIntLineEdit( m_Lines[i] );
    else if ( m_Lines[i].type == "DoubleLineEdit" ) addDoubleLineEdit( m_Lines[i] );
    else if ( m_Lines[i].type == "TextLineEdit" ) addTextLineEdit( m_Lines[i] );
    else if ( m_Lines[i].type == "CheckBox" ) addCheckBox( m_Lines[i] );
    else if ( m_Lines[i].type == "SpinBox" ) addSpinBox( m_Lines[i] );
    else if ( m_Lines[i].type == "DoubleSpinBox" ) addDoubleSpinBox( m_Lines[i] );
    else qDebug() << "Unknown type";
  }
}

void TemplateFormLayout::addComboBox( TemplateLine line )
{
  qDebug() << "Adding a ComboBox...";
  QComboBox* combobox = new QComboBox;
  QStringList description;
  QStringList value;
  QStringList L_open = line.options.split( "(" );
  for ( int i = 1; i < L_open.size(); i++ ) {
    QStringList L_close = L_open[i].split( ")" );
    QStringList L_elements = L_close[0].split( "," );
    description << L_elements[0];
    value << L_elements[1].trimmed();
  }
  int current = value.indexOf(line.getDefaultValue().trimmed());
  combobox->addItems( description );
  combobox-> setCurrentIndex(current);
  this->addRow( line.name, combobox );
  m_ComboBoxVector.push_back( combobox );
  m_ComboboxValues.push_back( value );
}

void TemplateFormLayout::addIntLineEdit( TemplateLine line )
{
  qDebug() << "Adding a IntLineEdit...";
  QValidator *validator = new QIntValidator( this );
  QLineEdit* int_lineedit = new QLineEdit;
  int_lineedit->setValidator( validator );
  line.default_value_of = line.options;
  int_lineedit->setText( line.getDefaultValue().trimmed() );
  this->addRow( line.name, int_lineedit );
  m_IntLineEditVector.push_back( int_lineedit );
}

void TemplateFormLayout::addDoubleLineEdit( TemplateLine line )
{
  qDebug() << "Adding a DoubleLineEdit...";
  QValidator *validator = new QDoubleValidator( this );
  QLineEdit* double_lineedit = new QLineEdit;
  double_lineedit->setValidator( validator );
  line.default_value_of = line.options;
  double_lineedit->setText( line.getDefaultValue().trimmed() );
  this->addRow( line.name, double_lineedit );
  m_DoubleLineEditVector.push_back( double_lineedit );
}

void TemplateFormLayout::addTextLineEdit( TemplateLine line )
{
  qDebug() << "Adding a TextLineEdit...";
  QLineEdit* text_lineedit = new QLineEdit;
  line.default_value_of = line.options;
  text_lineedit->setText( line.getDefaultValue().trimmed() );
  this->addRow( line.name, text_lineedit );
  m_TextLineEditVector.push_back( text_lineedit );
}

void TemplateFormLayout::addCheckBox( TemplateLine line )
{
  qDebug() << "Adding a CheckBox...";
  QCheckBox* check_box = new QCheckBox;
  QStringList L = line.options.split( "," );
  line.default_value_of = L[0];
  if ( line.getDefaultValue().trimmed() == "checked" ) check_box->setCheckState( Qt::Checked );
  else check_box->setCheckState( Qt::Unchecked );
  QPair < QString, QString > values;
  values.first = L[1];
  values.second = L[2];
  this->addRow( line.name, check_box );
  m_CheckBoxVector.push_back( check_box );
  m_CheckBoxValues.push_back( values );
}

void TemplateFormLayout::addSpinBox( TemplateLine line )
{
  qDebug() << "Adding a SpinBox...";
  QSpinBox* spin_box = new QSpinBox;
  QStringList L = line.options.split( "," );
  int minimum = L[0].trimmed().toInt();
  int maximum = L[1].trimmed().toInt();
  int step = L[2].trimmed().toInt();
  line.default_value_of = L[3];
  int value = line.getDefaultValue().trimmed().toInt();
  spin_box->setRange( minimum, maximum );
  spin_box->setSingleStep( step );
  spin_box->setValue( value );
  this->addRow( line.name, spin_box );
  m_SpinBoxVector.push_back( spin_box );
}

void TemplateFormLayout::addDoubleSpinBox( TemplateLine line )
{
  qDebug() << "Adding a DoubleSpinBox...";
  QDoubleSpinBox* double_spin_box = new QDoubleSpinBox;
  QStringList L = line.options.split( "," );
  double minimum = L[0].trimmed().toDouble();
  double maximum = L[1].trimmed().toDouble();
  double step = L[2].trimmed().toDouble();
  int decimals = L[3].trimmed().toInt();
  line.default_value_of = L[4];
  double value = line.getDefaultValue().trimmed().toDouble();
  double_spin_box->setRange( minimum, maximum );
  double_spin_box->setSingleStep( step );
  double_spin_box->setDecimals( decimals );
  double_spin_box->setValue( value );
  this->addRow( line.name, double_spin_box );
  m_DoubleSpinBoxVector.push_back( double_spin_box );
}

void TemplateFormLayout::open()
{
}

void TemplateFormLayout::save()
{
  qDebug() << "Saving...";
  getValues();
  m_file_template.save_egc();
}

void TemplateFormLayout::saveAs()
{
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
    if ( m_Lines[i].type == "ComboBox" ) {
      m_Lines[i].default_value_egc = readComboBox( combobox_idx );
      combobox_idx++;
    }
    else if ( m_Lines[i].type == "IntLineEdit" ) {
      m_Lines[i].default_value_egc =  readIntLineEdit( intlineedit_idx );
      intlineedit_idx++;
    }
    else if ( m_Lines[i].type == "DoubleLineEdit" ) {
      m_Lines[i].default_value_egc =  readDoubleLineEdit( doublelineedit_idx );
      doublelineedit_idx++;
    }
    else if ( m_Lines[i].type == "TextLineEdit" ) {
      m_Lines[i].default_value_egc =  readTextLineEdit( textlineedit_idx );
      textlineedit_idx++;
    }
    else if ( m_Lines[i].type == "CheckBox" ) {
      m_Lines[i].default_value_egc =  readCheckBox( checkbox_idx );
      checkbox_idx++;
    }
    else if ( m_Lines[i].type == "SpinBox" ) {
      m_Lines[i].default_value_egc =  readSpinBox( spinbox_idx );
      spinbox_idx++;
    }
    else if ( m_Lines[i].type == "DoubleSpinBox" ) {
      m_Lines[i].default_value_egc =  readDoubleSpinBox( doublespinbox_idx );
      doublespinbox_idx++;
    }
    else qDebug() << "Unknown type";
  }
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
