#include "dialoglineedit.h"

#include <QtGui>
#include <QHBoxLayout>
#include <QtDebug>
#include <QFileInfo>

#include "utilities.h"

DialogLineEdit::DialogLineEdit(QWidget *parent)
    : QWidget(parent)
{
  build();
}

DialogLineEdit::DialogLineEdit(bool openfile, QWidget *parent)
    : QWidget(parent)
{
  m_openfile = openfile;
  build();
}

void DialogLineEdit::build()
{
  setWindowTitle(tr("lineedit with file/directory selection dialog"));
//     this->setHeight(25);
  QHBoxLayout* layout = new QHBoxLayout();
  layout->setSpacing(-1);
  this->setContentsMargins(0, 0, 0, 0);
  layout->setContentsMargins(0, 0, 0, 0);
  this->setLayout(layout);
  m_button = new QPushButton(tr("..."));
  m_lineedit = new QLineEdit();
  m_button->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
  m_lineedit->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
  layout->addWidget(m_lineedit);
  layout->addWidget(m_button);

  connect(m_button, SIGNAL(clicked()), this, SLOT(openDialog()));

  m_button->setMaximumWidth(16);
  this->resize(180, 25);
//    this->resize(0, 0);
// setGeometry(10, 40, 180, 40);

  // by default, open directory and use file dialog
  m_openfile = false;
  m_UseFileDialog = true;
}

void DialogLineEdit::openDialog()
{
  emit clicked();

  if (m_UseFileDialog) {
    QString name;
    qDebug() << "m_dir=" << m_dir;
    if (m_openfile) name = QFileDialog::getOpenFileName(NULL, m_caption, m_dir, m_filter);
    else name = getDirectory(NULL, m_caption, m_dir);

    if (!name.isNull()) {
      m_lineedit->setText(name);
      m_dir = name;//getDir(name);
    }
  }
}

QString DialogLineEdit::text()
{
  return(m_lineedit->text());
}

void DialogLineEdit::setText(QString str)
{
  m_lineedit->setText(str);
  m_dir = str;//getDir(str);
}

bool DialogLineEdit::getMode()
{
  return m_openfile;
}

void DialogLineEdit::setMode(bool mode)
{
  m_openfile = mode;
}

void DialogLineEdit::setCaption(QString str)
{
  m_caption = str;
}

void DialogLineEdit::setDir(QString str)
{
  m_dir = str;
}

void DialogLineEdit::setFilter(QString str)
{
  m_filter = str;
}

QString DialogLineEdit::getCaption()
{
  return m_caption;
}

QString DialogLineEdit::getDir()
{
  return m_dir;
}

QString DialogLineEdit::getFilter()
{
  return m_filter;
}

// QString DialogLineEdit::getDir(QString str)
// {
//   if (m_openfile) {
//     return str;
//   }
//   else {
//     QFileInfo fileinfo(str);
//     return fileinfo.absolutePath();
//   }
// }

void DialogLineEdit::useFileDialog(bool arg)
{
  m_UseFileDialog = arg;
}
