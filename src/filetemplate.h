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
  
    void addComboBox(TemplateLine line);
    void getValues();
    QString readComboBox(int idx);
    
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
    FileTemplate file_template;
    QVector <TemplateLine> m_Lines;
    QStringList m_OutValues;
    QVector <QComboBox*> m_ComboBoxVector;
    QVector <QStringList> m_ComboboxValues;
};

#endif
