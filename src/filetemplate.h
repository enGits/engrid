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
#ifndef FILETEMPLATE_H
#define FILETEMPLATE_H

#include <QString>
#include <QVector>
#include <QFileInfo>
#include <QDialog>
#include <QtGui>
#include <QApplication>
#include <QStringList>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPair>

int fileTemplateTest( int argc, char ** argv );
int fileTemplateTest();

///@@@ TODO: Improve class structure

class TemplateLine
{
  public:
    QString m_Type;
    QString m_Name;
    QString m_Options;
    QString m_DefaultValueEgc;
    QString m_DefaultValueOpenFOAM;
    int m_Position;
  public:
    void print();
    QString getDefaultValue();
};

class FileTemplate
{
  public:
    QFileInfo m_FileInfo;
    QVector <TemplateLine> m_Lines;
    QString m_InText;
    QString m_OutText;
    QString m_Section;

  public:
    FileTemplate();
    FileTemplate( QString filename, QString section );
    int process();
    void print();
    QVector <TemplateLine> getLines();
    void setLines( QVector <TemplateLine> lines );
    QString getContents();
    void setContents( QString contents );
  public:
    int open( QString filename, QString section );
    int saveEgc();
    int exportToOpenFOAM( QString filename );
};

class TemplateFormLayout : public QFormLayout
{
    Q_OBJECT;

  private:
    QVector <FileTemplate> m_FileTemplate;
//     QVector < QVector <TemplateLine> > m_Lines;

    QVector <QComboBox*> m_ComboBoxVector;
    QVector <QStringList> m_ComboboxValues;
    QVector <QLineEdit*> m_IntLineEditVector;
    QVector <QLineEdit*> m_DoubleLineEditVector;
    QVector <QLineEdit*> m_TextLineEditVector;
    QVector <QCheckBox*> m_CheckBoxVector;
    QVector < QPair <QString, QString> > m_CheckBoxValues;
    QVector <QSpinBox*> m_SpinBoxVector;
    QVector <QDoubleSpinBox*> m_DoubleSpinBoxVector;

  public:
    TemplateFormLayout( QVector <QString> filename, QString section, char *name = 0, QWidget *parent = 0 );

    void addComboBox( TemplateLine line );
    void addIntLineEdit( TemplateLine line );
    void addDoubleLineEdit( TemplateLine line );
    void addTextLineEdit( TemplateLine line );
    void addCheckBox( TemplateLine line );
    void addSpinBox( TemplateLine line );
    void addDoubleSpinBox( TemplateLine line );

    QString readComboBox( int idx );
    QString readIntLineEdit( int idx );
    QString readDoubleLineEdit( int idx );
    QString readTextLineEdit( int idx );
    QString readCheckBox( int idx );
    QString readSpinBox( int idx );
    QString readDoubleSpinBox( int idx );

  public:
    void saveEgc();
};

class TemplateDialog : public QDialog
{
    Q_OBJECT;

  public:
    TemplateDialog( QVector <QString> files, QString section, QWidget *parent = 0 );
  private:
    QVector <QFileInfo> m_FileInfo;
    QVector <TemplateFormLayout*> m_TemplateFormLayoutVector;
  private slots:
    void saveEgc();
};

#endif
