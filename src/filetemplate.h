#ifndef FILETEMPLATE_H
#define FILETEMPLATE_H

#include <QString>
#include <QVector>
#include <QFileInfo>
#include <QDialog>
#include <QtGui>
#include <QApplication>
#include <QStringList>

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
  
public:
  FileTemplate();
  FileTemplate(QString filename);
  int open(QString filename);
  int save();
  int saveAs(QString filename);
  int process();
  void print();
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

private slots:
  void open();
  void save();
  void saveAs();
  
private:
  QPushButton *openButton;
  QPushButton *saveButton;
  QPushButton *saveAsButton;
};

#endif
