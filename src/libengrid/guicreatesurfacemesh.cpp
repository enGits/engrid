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

  connect(m_Ui.m_CheckBoxBoundaryLayerPolyhedral, SIGNAL(toggled(bool)), this, SLOT(updateTriangulateFarfieldState(bool)));
  updateTriangulateFarfieldState(m_Ui.m_CheckBoxBoundaryLayerPolyhedral->isChecked());

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

  connect(m_Ui.pushButton_SelectAll_BC,   SIGNAL(clicked()), this, SLOT(SelectAll_BC()));
  connect(m_Ui.pushButton_ClearAll_BC,    SIGNAL(clicked()), this, SLOT(ClearAll_BC()));
  connect(m_Ui.pushButton_SelectAll_BC_2, SIGNAL(clicked()), this, SLOT(SelectAll_BLayer()));
  connect(m_Ui.pushButton_ClearAll_BC_2,  SIGNAL(clicked()), this, SLOT(ClearAll_BLayer()));

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


void GuiCreateSurfaceMesh::updateTriangulateFarfieldState(bool on)
{
  m_Ui.m_CheckBoxTriangulateFarfield->setEnabled(on);
  m_Ui.m_CheckBoxMergeTopLayer->setEnabled(on);
  if (!on) {
    m_Ui.m_CheckBoxTriangulateFarfield->setChecked(false);
    m_Ui.m_CheckBoxMergeTopLayer->setChecked(false);
  }
}

///////////////////////////////////////////

int GuiCreateSurfaceMesh::readSettings()
{
  // Surface settings
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings");
  if (!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    QString first_token;
    in >> first_token;
    bool new_format = (first_token == "v2");

    QString max_edge, min_edge, growth;
    double nodes_per_quarter_circle;
    double feat2d, feat3d;

    if (new_format) {
      in >> max_edge >> min_edge >> growth >> nodes_per_quarter_circle;
    } else {
      // old format: first_token is max_edge
      max_edge = first_token;
      in >> min_edge >> growth >> nodes_per_quarter_circle;
    }
    m_Ui.lineEditMaximalEdgeLength->setText(max_edge);
    m_Ui.lineEditMinimalEdgeLength->setText(min_edge);
    m_Ui.lineEditGrowthFactor->setText(growth);
    m_Ui.doubleSpinBoxCurvature->setValue(nodes_per_quarter_circle);

    int num_bcs;
    in >> num_bcs;
    if (new_format) {
      // New format: each line is "bc_code state"
      QMap<int, int> bc_state_map;
      for (int i = 0; i < num_bcs; ++i) {
        int code, state;
        in >> code >> state;
        bc_state_map[code] = state;
      }
      for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
        int bc = m_Ui.listWidget->item(i)->text().split(" ")[0].toInt();
        if (bc_state_map.contains(bc)) {
          m_Ui.listWidget->item(i)->setCheckState(bc_state_map[bc] ? Qt::Checked : Qt::Unchecked);
        }
      }
    } else {
      // Old format: states in order, only restore if count matches
      if (num_bcs == m_Ui.listWidget->count()) {
        for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
          int check_state;
          in >> check_state;
          m_Ui.listWidget->item(i)->setCheckState(check_state ? Qt::Checked : Qt::Unchecked);
        }
      }
    }
    in >> feat2d;
    in >> feat3d;
    m_Ui.doubleSpinBox2DFeature->setValue(feat2d);
    m_Ui.m_DoubleSpinBox3DFeature->setValue(feat3d);
  }

  m_ELSManager.read();

  // Boundary layer settings
  buffer = GuiMainWindow::pointer()->getXmlSection("engrid/blayer/settings");
  if (!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    QString first_token;
    in >> first_token;
    bool new_format = (first_token == "v2");

    int num_bcs;
    if (new_format) {
      in >> num_bcs;
    } else {
      num_bcs = first_token.toInt();
    }

    if (new_format) {
      // New format: each line is "bc_code state"
      QMap<int, int> bc_state_map;
      for (int i = 0; i < num_bcs; ++i) {
        int code, state;
        in >> code >> state;
        bc_state_map[code] = state;
      }
      for (int i = 0; i < m_Ui.m_ListWidgetPrismaticLayers->count(); ++i) {
        int bc = m_Ui.m_ListWidgetPrismaticLayers->item(i)->text().split(" ")[0].toInt();
        if (bc_state_map.contains(bc)) {
          m_Ui.m_ListWidgetPrismaticLayers->item(i)->setCheckState(bc_state_map[bc] ? Qt::Checked : Qt::Unchecked);
        }
      }
    } else {
      // Old format: states in order
      if (num_bcs == m_Ui.m_ListWidgetPrismaticLayers->count()) {
        for (int i = 0; i < m_Ui.m_ListWidgetPrismaticLayers->count(); ++i) {
          int check_state;
          in >> check_state;
          m_Ui.m_ListWidgetPrismaticLayers->item(i)->setCheckState(check_state ? Qt::Checked : Qt::Unchecked);
        }
      }
    }

    double dv;
    int iv;
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerFeatureAngle->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerStretchingRatio->setValue(dv);

    if (new_format) {
      // New params: first layer height and number of layers
      QString h0_str;
      in >> h0_str;
      m_Ui.m_LineEditFirstLayerHeight->setText(h0_str);
      in >> iv; m_Ui.m_SpinBoxNumLayers->setValue(iv);
    } else {
      // Old format had farfield ratio here — skip it, use defaults
      double farfield_ratio;
      in >> farfield_ratio;
      // Leave first layer height and num layers at UI defaults
    }

    in >> iv; m_Ui.m_SpinBoxNormalRelaxationIterations->setValue(iv);
    in >> iv; m_Ui.m_SpinBoxHeightRelaxationIterations->setValue(iv);

    { double passband_dummy; in >> passband_dummy; } // passband no longer used

    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerLowerFaceLimit->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerUpperFaceLimit->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerFaceAngle->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerGapHeight->setValue(dv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerRadarAngle->setValue(dv);
    in >> iv; m_Ui.m_CheckBoxBoundaryLayerNormalVectorGrouping->setChecked(iv);
    in >> dv; m_Ui.m_DoubleSpinBoxBoundaryLayerGroupingAngle->setValue(dv);
    in >> iv; m_Ui.m_SpinBoxBufferLayers->setValue(iv);
    if (!in.atEnd()) {
      in >> iv;
      m_Ui.m_CheckBoxBoundaryLayerPolyhedral->setChecked(iv);
    } else {
      m_Ui.m_CheckBoxBoundaryLayerPolyhedral->setChecked(false);
    }
    if (!in.atEnd()) {
      in >> iv;
      m_Ui.m_CheckBoxTriangulateFarfield->setChecked(iv);
    } else {
      m_Ui.m_CheckBoxTriangulateFarfield->setChecked(false);
    }
    if (!in.atEnd()) {
      in >> iv;
      m_Ui.m_CheckBoxMergeTopLayer->setChecked(iv);
    } else {
      m_Ui.m_CheckBoxMergeTopLayer->setChecked(false);
    }
    if (!in.atEnd()) {
      in >> dv;
      m_Ui.m_DoubleSpinBoxMaxShellAngle->setValue(dv);
    } else {
      m_Ui.m_DoubleSpinBoxMaxShellAngle->setValue(60.0);
    }
    if (!in.atEnd()) {
      in >> dv;
      m_Ui.m_DoubleSpinBoxSmoothShellA->setValue(dv);
    } else {
      m_Ui.m_DoubleSpinBoxSmoothShellA->setValue(0.5);
    }
    if (!in.atEnd()) {
      in >> dv;
      m_Ui.m_DoubleSpinBoxSmoothShellB->setValue(dv);
    } else {
      m_Ui.m_DoubleSpinBoxSmoothShellB->setValue(2.0);
    }
    if (!in.atEnd()) {
      in >> dv;
      m_Ui.m_SpinBoxSmoothShellIter->setValue(static_cast<int>(dv));
    } else {
      m_Ui.m_SpinBoxSmoothShellIter->setValue(50);
    }
    if (!in.atEnd()) {
      in >> dv;
      m_Ui.m_SpinBoxLaplaceShellSweeps->setValue(static_cast<int>(dv));
    } else {
      m_Ui.m_SpinBoxLaplaceShellSweeps->setValue(20);
    }
    if (!in.atEnd()) {
      in >> dv;
      m_Ui.m_DoubleSpinBoxLaplaceShellOmega->setValue(dv);
    } else {
      m_Ui.m_DoubleSpinBoxLaplaceShellOmega->setValue(0.5);
    }
    // Optional trailing fields are ignored if present
    while (!in.atEnd()) { double dummy; in >> dummy; }
  }

  return(0);
}

int GuiCreateSurfaceMesh::writeSettings()
{
  // Surface settings — new format v2 with BC codes
  QString buffer = "";
  {
    QTextStream out(&buffer, QIODevice::WriteOnly);
    out << "v2\n";
    out << m_Ui.lineEditMaximalEdgeLength->text() << "\n";
    out << m_Ui.lineEditMinimalEdgeLength->text() << "\n";
    out << m_Ui.lineEditGrowthFactor->text() << "\n";
    out << m_Ui.doubleSpinBoxCurvature->value() << "\n";
    out << m_Ui.listWidget->count() << "\n";
    for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
      int bc = m_Ui.listWidget->item(i)->text().split(" ")[0].toInt();
      int state = (m_Ui.listWidget->item(i)->checkState() == Qt::Checked) ? 1 : 0;
      out << bc << " " << state << "\n";
    }
    out << m_Ui.doubleSpinBox2DFeature->value() << "\n";
    out << m_Ui.m_DoubleSpinBox3DFeature->value() << "\n";
  }
  GuiMainWindow::pointer()->setXmlSection("engrid/surface/settings", buffer);

  // Boundary layer settings — new format v2
  buffer = "";
  {
    QTextStream out(&buffer, QIODevice::WriteOnly);
    out << "v2\n";
    out << m_Ui.m_ListWidgetPrismaticLayers->count() << "\n";
    for (int i = 0; i < m_Ui.m_ListWidgetPrismaticLayers->count(); ++i) {
      int bc = m_Ui.m_ListWidgetPrismaticLayers->item(i)->text().split(" ")[0].toInt();
      int state = (m_Ui.m_ListWidgetPrismaticLayers->item(i)->checkState() == Qt::Checked) ? 1 : 0;
      out << bc << " " << state << "\n";
    }
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerFeatureAngle->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerStretchingRatio->value() << "\n";
    out << m_Ui.m_LineEditFirstLayerHeight->text() << "\n";
    out << m_Ui.m_SpinBoxNumLayers->value() << "\n";
    out << m_Ui.m_SpinBoxNormalRelaxationIterations->value() << "\n";
    out << m_Ui.m_SpinBoxHeightRelaxationIterations->value() << "\n";
    out << 0.0 << "\n"; // passband placeholder (field kept for file format compatibility)
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerLowerFaceLimit->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerUpperFaceLimit->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerFaceAngle->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerGapHeight->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerRadarAngle->value() << "\n";
    out << (m_Ui.m_CheckBoxBoundaryLayerNormalVectorGrouping->checkState() == Qt::Checked ? 1 : 0) << "\n";
    out << m_Ui.m_DoubleSpinBoxBoundaryLayerGroupingAngle->value() << "\n";
    out << m_Ui.m_SpinBoxBufferLayers->value() << "\n";
    out << (m_Ui.m_CheckBoxBoundaryLayerPolyhedral->checkState() == Qt::Checked ? 1 : 0) << "\n";
    out << (m_Ui.m_CheckBoxTriangulateFarfield->checkState() == Qt::Checked ? 1 : 0) << "\n";
    out << (m_Ui.m_CheckBoxMergeTopLayer->checkState() == Qt::Checked ? 1 : 0) << "\n";
    out << m_Ui.m_DoubleSpinBoxMaxShellAngle->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxSmoothShellA->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxSmoothShellB->value() << "\n";
    out << m_Ui.m_SpinBoxSmoothShellIter->value() << "\n";
    out << m_Ui.m_SpinBoxLaplaceShellSweeps->value() << "\n";
    out << m_Ui.m_DoubleSpinBoxLaplaceShellOmega->value() << "\n";
  }
  GuiMainWindow::pointer()->setXmlSection("engrid/blayer/settings", buffer);

  m_ELSManager.write();

  GuiMainWindow::pointer()->setXmlSection("engrid/surface/rules", m_Ui.textEdit->toPlainText());

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

void GuiCreateSurfaceMesh::SelectAll_BLayer()
{
  for (int i = 0; i < m_Ui.m_ListWidgetPrismaticLayers->count(); ++i) {
    m_Ui.m_ListWidgetPrismaticLayers->item(i)->setCheckState(Qt::Checked);
  }
}

void GuiCreateSurfaceMesh::ClearAll_BLayer()
{
  for (int i = 0; i < m_Ui.m_ListWidgetPrismaticLayers->count(); ++i) {
    m_Ui.m_ListWidgetPrismaticLayers->item(i)->setCheckState(Qt::Unchecked);
  }
}

void GuiCreateSurfaceMesh::setTextFromTable()
{
  m_Ui.textEdit->setText(GuiMainWindow::pointer()->getXmlSection("engrid/surface/rules").trimmed());
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
}
