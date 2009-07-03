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

class TemplateLine
{
  public:
    QString type;
    QString name;
    QString options;
    QString default_value;
    int position;
  public:
    void print();
};

class FileTemplate
{
  public:
    QFileInfo m_FileInfo;
    QVector <TemplateLine> m_Lines;
    QString m_InText;
    QString m_OutText;
    QStringList m_OutValues;

  public:
    FileTemplate();
    FileTemplate( QString filename );
    int open( QString filename );
    int save();
    int saveAs( QString filename );
    int process();
    void print();
    QVector <TemplateLine> getLines();
    void setOutValues( QStringList L );
    QString getContents();
    void setContents(QString contents);
};

class TemplateFormLayout : public QFormLayout
{
    Q_OBJECT

  private:
    FileTemplate m_file_template;
    QVector <TemplateLine> m_Lines;
    QStringList m_OutValues;

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
    TemplateFormLayout( QString filename, char *name = 0, QWidget *parent = 0 );

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

    void getValues();

    void open();
    void save();
    void saveAs();
};

class TemplateDialog : public QDialog
{
    Q_OBJECT
  public:
    TemplateDialog( QVector <QString> files, QWidget *parent = 0 );
  private:
    QVector <QFileInfo> m_FileInfo;
    QVector <TemplateFormLayout*> m_TemplateFormLayoutVector;
  private slots:
    void open();
    void save();
    void saveAs();
};

#endif
