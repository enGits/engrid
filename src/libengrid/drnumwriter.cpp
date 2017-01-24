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

#include "drnumwriter.h"
#include "guimainwindow.h"
#include "vtkEgNormalExtrusion.h"
#include "padsurface.h"

DrNumWriter::DrNumWriter()
{
  m_BackupGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
}

QList<BoundaryCondition> DrNumWriter::getBcsOfType(QString type)
{
  QList<BoundaryCondition> bcs_of_type;
  QSet<int> bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
  foreach (int bc_code, bcs) {
    BoundaryCondition bc = GuiMainWindow::pointer()->getBC(bc_code);
    PhysicalBoundaryCondition pbc = GuiMainWindow::pointer()->getPhysicalBoundaryCondition(bc.getType());
    if (pbc.getType() == type) {
      bcs_of_type << bc;
    }
  }
  qSort(bcs_of_type);
  return bcs_of_type;
}

void DrNumWriter::readSettings()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings");
  if(!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    in >> m_MaximalEdgeLength;
    in >> m_MinimalEdgeLength;
    in >> m_GrowthFactor;
  }
}

double DrNumWriter::edgeLength(QString bc_name)
{
  QString rules_txt = GuiMainWindow::pointer()->getXmlSection("engrid/surface/rules");
  rules_txt = rules_txt.replace("\n", " ");
  rules_txt = rules_txt.trimmed();
  QStringList rules = rules_txt.split(";", QString::SkipEmptyParts);
  double h_min = EG_LARGE_REAL;
  foreach (QString rule, rules) {
    QStringList parts = rule.split("=");
    if (parts.size() > 1) {
      QString left  = parts[0].trimmed();
      double h = parts[1].trimmed().toDouble();
      QStringList or_parts = left.split("<OR>");
      foreach (QString or_part, or_parts) {
        or_part = or_part.trimmed();
        QStringList and_parts = or_part.split("<AND>");
        if (and_parts.size() == 1) {
          QString and_part = and_parts[0].trimmed();
          if (and_part == bc_name) {
            h_min = min(h, h_min);
          }
        }
      }
    }
  }
  return h_min;
}

void DrNumWriter::prepareLevelSets(QList<BoundaryCondition> bcs, double distance)
{
  bool extrude = distance > 1e-3;
  foreach (BoundaryCondition bc, bcs) {
    BoundaryCondition new_bc = GuiMainWindow::pointer()->getBC(BoundaryCondition(bc.getName(), "level-set"));
    BoundaryCondition tmp_bc = GuiMainWindow::pointer()->getBC(BoundaryCondition(bc.getName(), "auxilary"));
    MeshPartition part(m_Grid);
    part.setBC(bc.getCode());
    if (!part.isPlanar()) {
      PhysicalBoundaryCondition pbc = GuiMainWindow::pointer()->getPhysicalBoundaryCondition(bc.getType());
      QString msg =  "only planar surfaces allowed for boundaries of type \"" + pbc.getType() + "\"";
      EG_ERR_RETURN(msg);
    }
    double A, P, Dh;
    vec3_t x, n;
    part.calcPlanarSurfaceMetrics(Dh, A, P, x, n);
    if (extrude) {
      QList<double> h;
      h << distance*Dh + 10*m_MaximalEdgeLength;
      x += h.last()*n;
      part.extrude(n, h, bc, false, false, true, BoundaryCondition(), BoundaryCondition(), tmp_bc);
      part.setBC(tmp_bc.getCode());
    }
    part.duplicate();
    part.resetBC(new_bc.getName(), new_bc.getType());

    //part.setBC(new_bc.getCode());
    //part.scale((Dh + 10*m_MaximalEdgeLength)/Dh, x);

    PadSurface pad;
    pad.setGrid(m_Grid);
    pad.addBC(new_bc.getCode());
    pad.relativeOff();
    pad.setDistance(10*m_MaximalEdgeLength);
    pad.setNewBC(new_bc.getCode());
    pad();
    part.setBC(new_bc.getCode());

    if (extrude) {
      part.translate(-10*m_MaximalEdgeLength*n);
    }
    part.writeSTL(getFileName() + "/engrid/" + new_bc.getName() + ".stl");
    QFile info_file(getFileName() + "/engrid/" + new_bc.getName() + ".dnc");
    info_file.open(QIODevice::WriteOnly);
    QTextStream info(&info_file);
    PhysicalBoundaryCondition pbc = GuiMainWindow::pointer()->getPhysicalBoundaryCondition(bc.getType());
    int max_size = 6;
    for (int i = 0; i < pbc.getNumVars(); ++i) {
      max_size = max(max_size, pbc.getVarName(i).size());
    }
    info << "string " + QString("name").leftJustified(max_size) << " = " << bc.getName() << ";\n";
    info << "string " + QString("type").leftJustified(max_size) << " = " << pbc.getType() << ";\n";
    info << "real   " + QString("A").leftJustified(max_size) << " = " << A << ";\n";
    info << "real   " + QString("P").leftJustified(max_size) << " = " << P << ";\n";
    info << "real   " + QString("Dh").leftJustified(max_size) << " = " << Dh << ";\n";
    info << "vec3_t " + QString("centre").leftJustified(max_size) << " = (" << x[0] << ", " << x[1] << ", " << x[2] << ");\n";
    info << "vec3_t " + QString("normal").leftJustified(max_size) << " = (" << n[0] << ", " << n[1] << ", " << n[2] << ");\n";
    for (int i = 0; i < pbc.getNumVars(); ++i) {
      info << pbc.getVarType(i).leftJustified(7) << pbc.getVarName(i).leftJustified(max_size) << " = " << pbc.getVarValueAsString(i) << ";\n";
    }
    double h = edgeLength(bc.getName());
    if (h < m_MaximalEdgeLength) {
      info << "real   " + QString("h").leftJustified(max_size) << " = " << h << ";\n";
    }
  }
}

void DrNumWriter::prepareWallLevelSets(QList<BoundaryCondition> bcs)
{
  foreach (BoundaryCondition bc, bcs) {
    MeshPartition part(m_Grid);
    part.setBC(bc.getCode());
    part.writeSTL(getFileName() + "/engrid/" + bc.getName() + ".stl");
    QFile info_file(getFileName() + "/engrid/" + bc.getName() + ".dnc");
    info_file.open(QIODevice::WriteOnly);
    QTextStream info(&info_file);
    PhysicalBoundaryCondition pbc = GuiMainWindow::pointer()->getPhysicalBoundaryCondition(bc.getType());
    int max_size = 4;
    for (int i = 0; i < pbc.getNumVars(); ++i) {
      max_size = max(max_size, pbc.getVarName(i).size());
    }
    info << "string " + QString("name").leftJustified(max_size) << " = " << bc.getName() << ";\n";
    info << "string " + QString("type").leftJustified(max_size) << " = " << pbc.getType() << ";\n";
    for (int i = 0; i < pbc.getNumVars(); ++i) {
      info << pbc.getVarType(i).leftJustified(7) << pbc.getVarName(i).leftJustified(max_size) << " = " << pbc.getVarValueAsString(i) << ";\n";
    }
    double h = edgeLength(bc.getName());
    if (h < m_MaximalEdgeLength) {
      info << "real   " + QString("h").leftJustified(max_size) << " = " << h << ";\n";
    }
  }
}

void DrNumWriter::computeBBox()
{
  m_X1 = vec3_t(EG_LARGE_REAL, EG_LARGE_REAL, EG_LARGE_REAL);
  m_X2 = -1*m_X1;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    for (int i = 0; i < 3; ++i) {
      m_X1[i] = min(m_X1[i], x[i]);
      m_X2[i] = max(m_X2[i], x[i]);
    }
  }
}

void DrNumWriter::writeGlobals()
{
  QFile info_file(getFileName() + "/engrid/global.dnc");
  info_file.open(QIODevice::WriteOnly);
  QTextStream info(&info_file);

  info << "string name  = global;\n";
  info << "real   h_max = " << m_MaximalEdgeLength << ";\n";
  info << "real   gf    = " << m_GrowthFactor << ";\n";
  info << "vector x1    = (" << m_X1[0] << ", " << m_X1[1] << ", " << m_X1[2] << ");\n";
  info << "vector x2    = (" << m_X2[0] << ", " << m_X2[1] << ", " << m_X2[2] << ");\n";
}

void DrNumWriter::backup()
{
  makeCopy(m_Grid, m_BackupGrid);
}

void DrNumWriter::restore()
{
  makeCopy(m_BackupGrid, m_Grid);
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
  updateNodeIndex(m_Grid);
  updateCellIndex(m_Grid);
}

void DrNumWriter::operate()
{
  readOutputDirectory();  
  if (isValid()) {
    backup();
    readSettings();

    QList<BoundaryCondition> turb_duct_in_bcs  = getBcsOfType("turbulent-duct-inlet");
    QList<BoundaryCondition> lam_duct_in_bcs   = getBcsOfType("laminar-duct-inlet");
    QList<BoundaryCondition> lam_ext_in_bcs    = getBcsOfType("inflow");
    QList<BoundaryCondition> out_bcs           = getBcsOfType("outflow");
    QList<BoundaryCondition> cyclic_bcs        = getBcsOfType("cyclic");
    QList<BoundaryCondition> sym_bcs           = getBcsOfType("symmetry");
    QList<BoundaryCondition> turb_wall_bcs     = getBcsOfType("DrNUM-turbulent-wall");
    QList<BoundaryCondition> lam_wall_bcs      = getBcsOfType("wall");
    QList<BoundaryCondition> slip_wall_bcs     = getBcsOfType("inviscid-wall");

    QString root_path = getFileName();
    QDir    root_dir(root_path);
    QString engrid_path = "engrid";
    QDir    engrid_dir(root_path + "/" + engrid_path);
    if (!engrid_dir.exists()) {
      root_dir.mkdir(engrid_path);
    } else {
      // delete all files in the levelset directory
      QStringList files = engrid_dir.entryList(QDir::Files);
      foreach (QString file, files) {
        engrid_dir.remove(file);
      }
    }

    prepareLevelSets(turb_duct_in_bcs,  3.0);
    prepareLevelSets(lam_duct_in_bcs,   0.0);
    prepareLevelSets(lam_ext_in_bcs,    0.0);
    prepareLevelSets(out_bcs,           0.0);
    prepareLevelSets(cyclic_bcs,        0.0);
    prepareLevelSets(sym_bcs,           0.0);

    prepareWallLevelSets(turb_wall_bcs);
    prepareWallLevelSets(lam_wall_bcs);
    prepareWallLevelSets(slip_wall_bcs);

    computeBBox();
    writeGlobals();
    writeSolverParameters(getFileName());
    restore();
  }
}
