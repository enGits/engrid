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

class TemplateLine {
public:
  QString type;
  QString name;
  QString options;
  int position;
public:
  void print();
};

class FileTemplate{
public:
  QFileInfo m_FileInfo;
  QVector <TemplateLine> m_Lines;
  QString m_InText;
  QString m_OutText;
  QStringList m_OutValues;
  
public:
  FileTemplate();
  FileTemplate(QString filename);
  int open(QString filename);
  int save();
  int saveAs(QString filename);
  int process();
  void print();
  QVector <TemplateLine> getLines();
  void setOutValues(QStringList L);
};

class MultipleFileTemplate{
public:
  QVector <QFileInfo> m_FileInfos;
public:
  void addFile(QString filename);
};

class GuiTemplateViewer : public QDialog
{
  Q_OBJECT
  public:
    // constructors
    /**
    * Constructor
    * @param parent Parent QWidget
    */
    GuiTemplateViewer(QWidget *parent = 0);
    
    /**
    * Constructor
    * @param filename name of the template file to use
    * @param parent Parent QWidget
    */
    GuiTemplateViewer(QString filename, QWidget *parent = 0);
  
    /**
    * Constructor
    * @param multiple_file_template MultipleFileTemplate object to use
    * @param parent Parent QWidget
    */
    GuiTemplateViewer(MultipleFileTemplate multiple_file_template, QWidget *parent = 0);
  
    void addComboBox(TemplateLine line);
    void addIntLineEdit(TemplateLine line);
    void addDoubleLineEdit(TemplateLine line);
    void addTextLineEdit(TemplateLine line);
    void addCheckBox(TemplateLine line);
    void addSpinBox(TemplateLine line);
    void addDoubleSpinBox(TemplateLine line);
  
    void getValues();
    void addFile(QString filename);
  
    QString readComboBox(int idx);
    QString readIntLineEdit(int idx);
    QString readDoubleLineEdit(int idx);
    QString readTextLineEdit(int idx);
    QString readCheckBox(int idx);
    QString readSpinBox(int idx);
    QString readDoubleSpinBox(int idx);
  
  private slots:
    void open();
    void save();
    void saveAs();
  
  private:
    QPushButton *openButton;
    QPushButton *saveButton;
    QPushButton *saveAsButton;
    QFormLayout *formLayout;
    QVBoxLayout *mainLayout;
    QHBoxLayout *bottomLayout;
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
  
    MultipleFileTemplate m_Files;
};

class SuperBox : public QFormLayout
{
  Q_OBJECT
    private:
    FileTemplate m_file_template;
    QVector <TemplateLine> m_Lines;
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
    SuperBox(QString filename, char *name=0, QWidget *parent=0);
  
  void addComboBox(TemplateLine line);
  void addIntLineEdit(TemplateLine line);
  void addDoubleLineEdit(TemplateLine line);
  void addTextLineEdit(TemplateLine line);
  void addCheckBox(TemplateLine line);
  void addSpinBox(TemplateLine line);
  void addDoubleSpinBox(TemplateLine line);
  
  ~SuperBox()
  {
    qDebug() << "Deleted: " << QFormLayout::objectName();
  }
};

class SuperGui : public QDialog
{
  Q_OBJECT
public:
  // constructors
  /**
  * Constructor
  * @param parent Parent QWidget
  */
  SuperGui(MultipleFileTemplate multiple_file_template, QWidget *parent = 0);
private:
  MultipleFileTemplate m_Files;
};

#endif
