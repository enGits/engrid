// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#ifndef DIALOGLINEEDIT_H
#define DIALOGLINEEDIT_H

#include <QTime>
#include <QWidget>
#include <QtDesigner/QDesignerExportWidget>
#include <QLineEdit>
#include <QPushButton>

/// A QWidget with a QLineEdit and a small QPushButton.
/** It allows easy file or directory selection. The clicked() signal can also be used to call up a custom dialog instead. */
class
#ifndef WIN32
    QDESIGNER_WIDGET_EXPORT
#endif
    DialogLineEdit : public QWidget
{
    Q_OBJECT
    Q_PROPERTY(bool m_openfile READ getMode WRITE setMode)
    Q_PROPERTY(bool m_UseFileDialog READ getUseFileDialog WRITE useFileDialog)

  private:
    QPushButton *m_button;///< The QPushButton
    QLineEdit *m_lineedit;///< The QLineEdit
    QString m_caption;///< The caption to use in the dialog
    QString m_dir;///< The default directory to use in the dialog
    QString m_filter;///< The filter to use in the dialog
    bool m_openfile;///< if true a getOpenFileName dialog is used, otherwise a getExistingDirectory dialog is used.
    bool m_UseFileDialog;///< if true a standard dialog is used, otherwise nothing happens (but the clicked() signal can be used to open a custom dialog)

  private:
    void build();///< builds the widget
    //QString getDir(QString str);///< returns the path of the selected file or the selected path

  public:
    DialogLineEdit(QWidget *parent = 0);///< Default constructor
    /** Constructor.
    * @param openfile true = use open file dialog, false = use open directory dialog
    * @param parent parent widget
    */
    DialogLineEdit(bool openfile, QWidget *parent = 0);
    QString text();///< Get the text of the QLineEdit
    void setText(QString str);///< set the text of the QLineEdit

  public:
    void setCaption(QString str);///< set the caption to use in the file dialog
    void setDir(QString str);///< set the default directory to use in the file dialog
    void setFilter(QString str);///< set the filter to use in the file dialog
    QString getCaption();///< get the caption used in the file dialog
    QString getDir();///< get the default directory used in the file dialog
    QString getFilter();///< get the filter used in the file dialog
    bool getMode();///< get the mode: true = use open file dialog, false = use open directory dialog
    void setMode(bool mode);///< set the mode: true = use open file dialog, false = use open directory dialog
    void useFileDialog(bool arg);///< set whether or not to use standard dialogs when the button is clicked
    bool getUseFileDialog() {return m_UseFileDialog;}///< get whether or not to use standard dialogs when the button is clicked

  private slots:
    void openDialog();///< emit the clicked signal + open standard dialog if m_UseFileDialog is true

  signals:
    void clicked();///< This signal is emitted when the QPushButton is clicked
};

#endif
