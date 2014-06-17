// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#include "guicreatesurfacemesh.h"

#include "swaptriangles.h"
#include "surfacemesher.h"
#include "vertexdelegate.h"
#include "laplacesmoother.h"
#include "guimainwindow.h"
#include "containertricks.h"
#include "updatedesiredmeshdensity.h"
#include "utilities.h"

#include <vtkSmoothPolyDataFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkLongArray.h>
#include <vtkEgNormalExtrusion.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkCellLocator.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

#include <QtGui>
#include <QTextStream>

#include <iostream>
#include <fstream>
#include <cmath>

#include <cstdio>

GuiCreateSurfaceMesh::GuiCreateSurfaceMesh()
{
  EG_TYPENAME;
  
  setQuickSave(true);

  populateBoundaryCodes(m_Ui.listWidget);
  populateBoundaryCodes(m_Ui.m_ListWidgetPrismaticLayers);
  m_Ui.lineEditMaximalEdgeLength->setText("1000");

  //Load settings
  readSettings();
  
  Nbc = m_Ui.listWidget-> count ();
  
  current_filename= GuiMainWindow::pointer()->getFilename();

  {
    QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/table");
    QTextStream in(&buffer, QIODevice::ReadOnly);

    m_NumRows = 0;
    m_NumCols = 0;
    in >> m_NumRows >> m_NumCols;

    if (m_NumRows == 0) {
      m_NumCols = Nbc + 3;
    }
    if(m_NumCols != Nbc + 3) {
      //m_NumCols = Nbc + 3;
    }

    int row, column;
    QString str;

    m_Table.clear();
    m_Table.fill(QVector<QString>(m_NumCols), m_NumRows);

    for (int i = 0; i < m_NumRows; ++i) {
      for (int j = 0; j < m_NumCols; ++j) {
        in >> row >> column >> str;
        if (str == "{{{empty}}}") {
          str = "";
        }
        m_Table[row][column] = str;
      }
    }
  }
  setTextFromTable();
  m_Ui.m_TextEditPrismaticLayers->setText(GuiMainWindow::pointer()->getXmlSection("engrid/blayer/rules"));
  
  connect(m_Ui.pushButton_SelectAll_BC, SIGNAL(clicked()), this, SLOT(SelectAll_BC()));
  connect(m_Ui.pushButton_ClearAll_BC, SIGNAL(clicked()), this, SLOT(ClearAll_BC()));

  m_ELSManager.setListWidget(m_Ui.listWidgetSources);
  m_ELSManager.read();
  m_ELSManager.populateListWidget();
  connect(m_Ui.pushButtonAddSphere,    SIGNAL(clicked()), this, SLOT(addSphere()));
  connect(m_Ui.pushButtonAddCone,      SIGNAL(clicked()), this, SLOT(addCone()));
  connect(m_Ui.pushButtonAddBox,       SIGNAL(clicked()), this, SLOT(addBox()));
  connect(m_Ui.pushButtonAddPipe,      SIGNAL(clicked()), this, SLOT(addPipe()));
  connect(m_Ui.pushButtonEditSource,   SIGNAL(clicked()), this, SLOT(edit()));
  connect(m_Ui.pushButtonDeleteSource, SIGNAL(clicked()), this, SLOT(remove()));

}

///////////////////////////////////////////

int GuiCreateSurfaceMesh::readSettings()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings");
  if(!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    QString str;
    in >> str;
    m_Ui.lineEditMaximalEdgeLength->setText(str);
    in >> str;
    m_Ui.lineEditMinimalEdgeLength->setText(str);
    in >> str;
    m_Ui.lineEditGrowthFactor->setText(str);
    double nodes_per_quarter_circle;
    in >> nodes_per_quarter_circle;
    m_Ui.doubleSpinBoxCurvature->setValue(nodes_per_quarter_circle);
    int num_bcs;
    in >> num_bcs;
    if (num_bcs == m_Ui.listWidget->count()) {
      int check_state;
      for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
        in >> check_state;
        if (check_state == 1) {
          m_Ui.listWidget->item(i)->setCheckState(Qt::Checked);
        } else {
          m_Ui.listWidget->item(i)->setCheckState(Qt::Unchecked);
        }
      }
    }
    double v;
    in >> v;
    m_Ui.doubleSpinBox2DFeature->setValue(v);
    in >> v;
    m_Ui.m_DoubleSpinBox3DFeature->setValue(v);
  }

  m_ELSManager.read();

  buffer = GuiMainWindow::pointer()->getXmlSection("engrid/blayer/settings");
  if(!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    int num_bcs;
    in >> num_bcs;
    if (num_bcs == m_Ui.m_ListWidgetPrismaticLayers->count()) {
      int check_state;
      for (int i = 0; i < m_Ui.m_ListWidgetPrismaticLayers->count(); ++i) {
        in >> check_state;
        if (check_state == 1) {
          m_Ui.m_ListWidgetPrismaticLayers->item(i)->setCheckState(Qt::Checked);
        } else {
          m_Ui.m_ListWidgetPrismaticLayers->item(i)->setCheckState(Qt::Unchecked);
        }
      }
    }
    double dv;
    int iv;
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerFeatureAngle->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerStretchingRatio->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerFarfieldRatio->setValue(dv);
    in >> iv; m_Ui.m_SpinBoxNormalRelaxationIterations->setValue(iv);
    in >> iv; m_Ui.m_SpinBoxHeightRelaxationIterations->setValue(iv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerLowerFaceLimit->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerUpperFaceLimit->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerFaceAngle->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerGapHeight->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerRadarAngle->setValue(dv);
    in >> iv; m_Ui.m_CheckBoxBoundaryLayerNormalVectorGrouping->setChecked(iv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerGroupingAngle->setValue(dv);
  }

  return(0);
}

int GuiCreateSurfaceMesh::writeSettings()
{
  QString buffer = "";
  {
    QTextStream out(&buffer, QIODevice::WriteOnly);
    out << "\n";
    out << m_Ui.lineEditMaximalEdgeLength->text() << "\n";
    out << m_Ui.lineEditMinimalEdgeLength->text() << "\n";
    out << m_Ui.lineEditGrowthFactor->text() << "\n";
    out << m_Ui.doubleSpinBoxCurvature->value() << "\n";
    out << m_Ui.listWidget->count() << "\n";
    for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
      if (m_Ui.listWidget->item(i)->checkState() == Qt::Checked) {
        out << "1 \n";
      } else {
        out << "0 \n";
      }
    }
    out << m_Ui.doubleSpinBox2DFeature->value() << "\n";
    out << m_Ui.m_DoubleSpinBox3DFeature->value() << "\n";
  }
  GuiMainWindow::pointer()->setXmlSection("engrid/surface/settings", buffer);

  buffer = "";
  {
    QTextStream out(&buffer, QIODevice::WriteOnly);
    out << "\n";
    out << m_Ui.m_ListWidgetPrismaticLayers->count() << "\n";
    for (int i = 0; i < m_Ui.m_ListWidgetPrismaticLayers->count(); ++i) {
      if (m_Ui.m_ListWidgetPrismaticLayers->item(i)->checkState() == Qt::Checked) {
        out << "1 \n";
      } else {
        out << "0 \n";
      }
    }
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerFeatureAngle->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerStretchingRatio->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerFarfieldRatio->value() << "\n";
    out << m_Ui.m_SpinBoxNormalRelaxationIterations->value() << "\n";
    out << m_Ui.m_SpinBoxHeightRelaxationIterations->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerLowerFaceLimit->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerUpperFaceLimit->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerFaceAngle->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerGapHeight->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerRadarAngle->value() << "\n";
    if (m_Ui.m_CheckBoxBoundaryLayerNormalVectorGrouping->checkState() == Qt::Checked) out << "1\n";
    else                                                                               out << "0\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerGroupingAngle->value() << "\n";
  }
  GuiMainWindow::pointer()->setXmlSection("engrid/blayer/settings", buffer);

  m_ELSManager.write();

  GuiMainWindow::pointer()->setXmlSection("engrid/surface/rules", m_Ui.textEdit->toPlainText());
  GuiMainWindow::pointer()->setXmlSection("engrid/blayer/rules", m_Ui.m_TextEditPrismaticLayers->toPlainText());

  return(0);
}

///////////////////////////////////////////

void GuiCreateSurfaceMesh::SelectAll_BC()
{
  for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
    m_Ui.listWidget->item(i)->setCheckState(Qt::Checked);
  }
}

void GuiCreateSurfaceMesh::ClearAll_BC()
{
  for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
    m_Ui.listWidget->item(i)->setCheckState(Qt::Unchecked);
  }
}

void GuiCreateSurfaceMesh::setTextFromTable()
{
  EG_STOPDATE("2014-07-01");
  m_Ui.textEdit->setText(GuiMainWindow::pointer()->getXmlSection("engrid/surface/rules"));
  m_Ui.m_TextEditPrismaticLayers->setText(GuiMainWindow::pointer()->getXmlSection("engrid/blayer/rules"));
  return;

  QString text = "";
  QSet<QString> values;
  for (int i = 0; i < m_NumRows; ++i) {
    values.insert(m_Table[i][m_NumCols-1]);
  }
  foreach (QString value, values) {
    int Ni = 0;
    for (int i = 0; i < m_NumRows; ++i) {
      if (m_Table[i][m_NumCols-1] == value) {
        ++Ni;
      }
    }
    int ni = 0;
    for (int i = 0; i < m_NumRows; ++i) {
      if (m_Table[i][m_NumCols-1] == value) {
        int Nj = 0;
        for (int j = 0; j < m_NumCols-3; ++j) {
          if (m_Table[i][j] == "2") {
            ++Nj;
          }
        }
        int nj = 0;
        for (int j = 0; j < m_NumCols-3; ++j) {
          if (m_Table[i][j] == "2") {
            QString bc = m_Ui.listWidget->item(j)->text().split(":")[1].trimmed();
            text += bc;
            ++nj;
            if (nj < Nj) {
              text += " <AND> ";
            }
          }
        }
        ++ni;
        if (ni < Ni) {
          text += " <OR>\n";
        }
      }
    }
    text += " = " + value + ";\n\n";
  }
  m_Ui.textEdit->setText(text);
}

void GuiCreateSurfaceMesh::getTableFromText()
{
  QMap<QString, int> bc_map;
  for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
    bc_map[m_Ui.listWidget->item(i)->text().split(":")[1].trimmed()] = i;
  }
  m_NumCols = bc_map.size() + 3;
  QStringList rules = m_Ui.textEdit->toPlainText().split(";", QString::SkipEmptyParts);
  m_Table.clear();
  foreach (QString rule, rules) {
    rule = rule.trimmed();
    QStringList parts = rule.split("=");
    if (parts.count() > 1) {
      QString left = parts[0].trimmed();
      QString right = parts[1].trimmed();
      QStringList rows = left.split("<OR>");
      foreach (QString row, rows) {
        row = row.trimmed();
        m_Table.append(QVector<QString>(m_NumCols, "1"));
        int r = m_Table.count() - 1;
        m_Table[r][m_NumCols - 3] = "any";
        m_Table[r][m_NumCols - 2] = "";
        m_Table[r][m_NumCols - 1] = right;
        QStringList cols = row.split("<AND>");
        foreach (QString col, cols) {
          col = col.trimmed();
          m_Table[r][bc_map[col]] = "2";
        }
      }
    }
  }
  m_NumRows = m_Table.size();
}

void GuiCreateSurfaceMesh::operate()
{
  writeSettings();
  getTableFromText();
  QString buffer = "";
  {
    QTextStream out(&buffer, QIODevice::WriteOnly);
    out << "\n";
    out << m_NumRows << " " << m_NumCols << "\n";
    for (int row = 0; row < m_NumRows; ++row) {
      for (int column = 0; column < m_NumCols; ++column) {
        QString str = m_Table[row][column];
        if (str.isEmpty()) {
          str = "{{{empty}}}";
        }
        out << row << " " << column << " " << str << "\n";
      }
    }
  }
  GuiMainWindow::pointer()->setXmlSection("engrid/surface/table", buffer);
  GuiMainWindow::pointer()->setXmlSection("engrid/blayer/rules", m_Ui.m_TextEditPrismaticLayers->toPlainText());
}
