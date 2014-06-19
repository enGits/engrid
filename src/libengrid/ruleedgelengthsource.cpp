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
// 
#include "ruleedgelengthsource.h"
#include "guimainwindow.h"

RuleEdgeLengthSource::RuleEdgeLengthSource(QString rule, vtkUnstructuredGrid *grid)
{
  readGrowthFactor();
  bool is_valid = true;
  QList<QList<int> > bc_combinations;
  rule = rule.trimmed();
  QList<QString> all_symbolic_bcs;
  QMap<QString, int> name2bc;
  QVector<int> bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
  foreach (int bc, bcs) {
    QString name = GuiMainWindow::pointer()->getBC(bc).getName();
    all_symbolic_bcs.append(name);
    name2bc[name] = bc;
  }
  QStringList parts = rule.split("=");
  if (parts.count() > 1) {
    QString left  = parts[0].trimmed();
    m_EdgeLength = parts[1].trimmed().toDouble();
    QStringList or_parts = left.split("<OR>");
    foreach (QString or_part, or_parts) {
      QList<int> bc_combination;
      or_part = or_part.trimmed();
      QStringList and_parts = or_part.split("<AND>");
      foreach (QString and_part, and_parts) {
        and_part = and_part.trimmed();
        if (!all_symbolic_bcs.contains(and_part)) {
          QString msg = "unknown boundary name \"" + and_part + "\" in rule:\"" + rule + "\"";
          EG_ERR_RETURN(msg);
        }
        int bc = name2bc[and_part];
        bc_combination.append(bc);
      }
      bc_combinations.append(bc_combination);
    }
  } else {
    is_valid = true;
  }
  if (!is_valid) {
    QString msg = "invalid rule:\"" + rule + "\"";
    EG_ERR_RETURN(msg);
  }
  MeshPartition part(grid, true);
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    QList<int> bcs;
    for (int i = 0; i < part.n2bcGSize(id_node); ++i) {
      bcs.append(part.n2bcG(id_node, i));
    }
    bool use_node = false;
    foreach (QList<int> bc_combination, bc_combinations) {
      bool found = true;
      foreach (int bc, bc_combination) {
        if (!bcs.contains(bc)) {
          found = false;
          break;
        }
      }
      if (found) {
        use_node = true;
        break;
      }
    }
    if (use_node) {
      vec3_t x;
      grid->GetPoint(id_node, x.data());
      m_Points.append(x);
    }
  }
  m_PointFinder.setPoints(m_Points);
}

void RuleEdgeLengthSource::readGrowthFactor()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings");
  if(!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    QString str;
    in >> str;
    in >> str;
    in >> str;
    m_GrowthFactor = str.toDouble();
  } else {
    m_GrowthFactor = 2.0;
  }
}

double RuleEdgeLengthSource::edgeLength(vec3_t x)
{
  QVector<int> close_points;
  m_PointFinder.getClosePoints(x, close_points);
  if (close_points.size() == 0) {
    EG_BUG;
  }
  double d = EG_LARGE_REAL;
  foreach (int i, close_points) {
    vec3_t xp = m_Points[i];
    d = min(d, (x - xp).abs());
  }
  double n = logarithm(m_GrowthFactor, 1.0 - (1.0 - m_GrowthFactor)*d/m_EdgeLength) - 1.0;
  return m_EdgeLength*pow(m_GrowthFactor, n);
}

