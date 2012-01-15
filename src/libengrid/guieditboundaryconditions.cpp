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

#include "guieditboundaryconditions.h"

#include "guimainwindow.h"
#include "volumedefinition.h"
#include "filetemplate.h"
#include "physicalboundarycondition.h"
#include "multipagewidget.h"
#include "multipagewidgetpage.h"

#include <QVBoxLayout>
#include <QFileInfo>

GuiEditBoundaryConditions::GuiEditBoundaryConditions()
{
  connect(ui.pushButtonAdd, SIGNAL(clicked()), this, SLOT(addVol()));
  connect(ui.pushButtonDelete, SIGNAL(clicked()), this, SLOT(delVol()));
  connect(ui.pushButtonAddBoundaryType, SIGNAL(clicked()), this, SLOT(addBoundaryType()));
  connect(ui.pushButtonDeleteBoundaryType, SIGNAL(clicked()), this, SLOT(deleteBoundaryType()));
  connect(ui.listWidgetBoundaryType, SIGNAL(itemSelectionChanged()), this, SLOT(changePhysicalValues()));
  connect(ui.pushButton_AddProcess, SIGNAL(clicked()), this, SLOT(addProcess()));
  connect(ui.pushButton_RemoveProcess, SIGNAL(clicked()), this, SLOT(deleteProcess()));
  connect(ui.pushButton_ImportHostFile, SIGNAL(clicked()), this, SLOT(importHostFile()));
  connect(ui.pushButton_ExportHostFile, SIGNAL(clicked()), this, SLOT(exportHostFile()));
  
  setupSolvers();
  loadMpiParameters();

  m_BcMap = NULL;
  delegate = new GuiVolumeDelegate();
  delegate->setFirstCol(3);

  //set initial tab
  ui.tabWidget->setCurrentIndex(0);
}

GuiEditBoundaryConditions::~GuiEditBoundaryConditions()
{
  delete delegate;
}

void GuiEditBoundaryConditions::before()
{
  if (!m_BcMap) EG_BUG;
  resetOrientation(m_Grid);
  while (ui.T->rowCount()) ui.T->removeRow(0);
  foreach(int i, m_BoundaryCodes) {
    BoundaryCondition bc = (*m_BcMap)[i];
    ui.T->insertRow(ui.T->rowCount());
    int r = ui.T->rowCount() - 1;
    ui.T->setItem(r, 0, new QTableWidgetItem());
    ui.T->item(r, 0)->setFlags(ui.T->item(r, 0)->flags() & (~Qt::ItemIsSelectable));
    ui.T->item(r, 0)->setFlags(ui.T->item(r, 0)->flags() & (~Qt::ItemIsEditable));
    ui.T->setItem(r, 1, new QTableWidgetItem());
    ui.T->setItem(r, 2, new QTableWidgetItem());
    QString idx;
    idx.setNum(i);
    ui.T->item(r, 0)->setText(idx);
    QString name = bc.getName();
    if (name == "unknown") name = QString("BC") + idx;
    ui.T->item(r, 1)->setText(name);
    ui.T->item(r, 2)->setText(bc.getType());
  }

  updateVol();
  updatePhysicalBoundaryConditions();

  ui.T->setItemDelegate(delegate);
}

void GuiEditBoundaryConditions::operate()
{
  // BoundaryCondition and VolumeDefinition
  QVector<VolumeDefinition> vols(ui.T->columnCount());
  for (int j = 3; j < ui.T->columnCount(); ++j) {
    QString vol_name = ui.T->horizontalHeaderItem(j)->text();
    VolumeDefinition V(vol_name, j - 2);
    vols[j] = V;
  }
  for (int i = 0; i < ui.T->rowCount(); ++i) {
    int bc = ui.T->item(i, 0)->text().toInt();
    BoundaryCondition BC(ui.T->item(i, 1)->text(), ui.T->item(i, 2)->text());
    (*m_BcMap)[bc] = BC;
    for (int j = 3; j < ui.T->columnCount(); ++j) {
      QString vol_name = ui.T->horizontalHeaderItem(j)->text();
      VolumeDefinition V = vols[j];
      if (ui.T->item(i, j)->text() == "green")  V.addBC(bc,  1);
      else if (ui.T->item(i, j)->text() == "yellow") V.addBC(bc, -1);
      else                                           V.addBC(bc,  0);
      vols[j] = V;
    }
  }
  QList<VolumeDefinition> vol_list;
  for (int j = 3; j < ui.T->columnCount(); ++j) {
    vol_list.append(vols[j]);
  }
  GuiMainWindow::pointer()->setAllVols(vol_list);

  // PhysicalBoundaryConditions
  savePhysicalValues();
  GuiMainWindow::pointer()->setAllPhysicalBoundaryConditions(m_PhysicalBoundaryConditionsMap);

  saveSolverParameters();
  saveMpiParameters();

}

//==========================

void GuiEditBoundaryConditions::updateVol()
{
  while (ui.T->columnCount() > 3) {
    ui.T->removeColumn(3);
  }
  QList<VolumeDefinition> vols = GuiMainWindow::pointer()->getAllVols();
  foreach(VolumeDefinition V, vols) {
    m_VolMap[V.getName()] = V;
    int c = ui.T->columnCount();
    ui.T->insertColumn(c);
    ui.T->setHorizontalHeaderItem(c, new QTableWidgetItem(V.getName()));
    for (int i = 0; i < ui.T->rowCount(); ++i) {
      int bc = ui.T->item(i, 0)->text().toInt();
      if (V.getSign(bc) ==  1) ui.T->setItem(i, c, new QTableWidgetItem("green"));
      else if (V.getSign(bc) == -1) ui.T->setItem(i, c, new QTableWidgetItem("yellow"));
      else                          ui.T->setItem(i, c, new QTableWidgetItem(" "));
    }
  }
}

void GuiEditBoundaryConditions::addVol()
{
  cout << "Adding volume" << endl;
  QString name = ui.lineEditVolume->text();
  if (!name.isEmpty() && !m_VolMap.contains(name)) {
    VolumeDefinition NV(name, ui.T->columnCount() - 2);
    m_VolMap[NV.getName()] = NV;
    int c = ui.T->columnCount();
    ui.T->insertColumn(c);
    ui.T->setHorizontalHeaderItem(c, new QTableWidgetItem(NV.getName()));
    for (int i = 0; i < ui.T->rowCount(); ++i) {
      int bc = ui.T->item(i, 0)->text().toInt();
      if (NV.getSign(bc) == 1)  ui.T->setItem(i, c, new QTableWidgetItem("green"));
      else if (NV.getSign(bc) == -1) ui.T->setItem(i, c, new QTableWidgetItem("yellow"));
      else                           ui.T->setItem(i, c, new QTableWidgetItem(" "));
    }
  }
}

void GuiEditBoundaryConditions::delVol()
{
  cout << "Deleting volume" << endl;
  int col = ui.T->currentColumn();
  cout << "col=" << col << endl;
  if (col > 2) {
    QString name = ui.T->horizontalHeaderItem(col)->text();
    cout << "name=" << qPrintable(name) << endl;
    m_VolMap.remove(name);
    ui.T->removeColumn(col);
  }
  else {
    cout << "Nothing to delete." << endl;
  }
}

//==========================

void GuiEditBoundaryConditions::loadPhysicalValues()
{
  while (ui.tableWidgetPBC->rowCount()) {
    ui.tableWidgetPBC->removeRow(0);
  }
  for (int i = 0; i < m_PBC_current.getNumVars(); ++i) {
    QString str;
    str.setNum(m_PBC_current.getVarValue(i));
    ui.tableWidgetPBC->insertRow(ui.tableWidgetPBC->rowCount());
    int r = ui.tableWidgetPBC->rowCount() - 1;
    ui.tableWidgetPBC->setItem(r, 0, new QTableWidgetItem());
    ui.tableWidgetPBC->item(r, 0)->setFlags(ui.tableWidgetPBC->item(r, 0)->flags() & (~Qt::ItemIsSelectable));
    ui.tableWidgetPBC->item(r, 0)->setFlags(ui.tableWidgetPBC->item(r, 0)->flags() & (~Qt::ItemIsEditable));
    ui.tableWidgetPBC->setItem(r, 1, new QTableWidgetItem());
    ui.tableWidgetPBC->item(r, 0)->setText(m_PBC_current.getVarName(i));
    ui.tableWidgetPBC->item(r, 1)->setText(str);
  }
  ui.tableWidgetPBC->resizeColumnsToContents();
}

void GuiEditBoundaryConditions::savePhysicalValues()
{
  if(m_PhysicalBoundaryConditionsMap.contains(m_PBC_current.getName())) {
    for (int i = 0; i < m_PBC_current.getNumVars(); ++i) {
      m_PBC_current.setValue(i, ui.tableWidgetPBC->item(i, 1)->text().toDouble());
    }
    m_PhysicalBoundaryConditionsMap[m_PBC_current.getName()] = m_PBC_current;
  }
}

void GuiEditBoundaryConditions::changePhysicalValues()
{
  if (ui.listWidgetBoundaryType->count() > 0) {
    int index = ui.listWidgetBoundaryType->currentRow();
    QString name = ui.listWidgetBoundaryType->currentItem()->text();
    savePhysicalValues();
    if(m_PhysicalBoundaryConditionsMap.contains(name)) {
      m_PBC_current = m_PhysicalBoundaryConditionsMap[name];
      m_PBC_current.setIndex(index);
    }
    else EG_BUG;
    loadPhysicalValues();
  }
}

void GuiEditBoundaryConditions::addBoundaryType()
{
  QString name = ui.lineEditBoundaryType->text();
  if (!name.isEmpty() && !m_PhysicalBoundaryConditionsMap.contains(name)) {
    PhysicalBoundaryCondition PBC;
    PBC.setName(ui.lineEditBoundaryType->text());
    PBC.setIndex(ui.listWidgetBoundaryType->count());
    PBC.setType(ui.comboBoxBoundaryType->currentText());
    m_PhysicalBoundaryConditionsMap[PBC.getName()] = PBC;
    ui.listWidgetBoundaryType->addItem(PBC.getName());
  }
}

void GuiEditBoundaryConditions::deleteBoundaryType()
{
  int row = ui.listWidgetBoundaryType->currentRow();
  if (ui.listWidgetBoundaryType->count() <= 1) {
    ui.listWidgetBoundaryType->clear();
    m_PhysicalBoundaryConditionsMap.clear();
  }
  else if (0 <= row && row < ui.listWidgetBoundaryType->count()) {
    QListWidgetItem* list_widget_item = ui.listWidgetBoundaryType->takeItem(row);
    m_PhysicalBoundaryConditionsMap.remove(list_widget_item->text());
    delete list_widget_item;
  }
}

void GuiEditBoundaryConditions::updatePhysicalBoundaryConditions()
{
  // clear list
  ui.listWidgetBoundaryType->clear();
  
  // fill list
  QList<PhysicalBoundaryCondition> physical_boundary_conditions = GuiMainWindow::pointer()->getAllPhysicalBoundaryConditions();
  foreach(PhysicalBoundaryCondition PBC, physical_boundary_conditions) {
    ui.listWidgetBoundaryType->addItem(PBC.getName());
    m_PhysicalBoundaryConditionsMap[PBC.getName()] = PBC;
  }
  
  // select and load first item
  if (ui.listWidgetBoundaryType->count() > 0) {
    ui.listWidgetBoundaryType->setCurrentRow(0);
    QString name = ui.listWidgetBoundaryType->currentItem()->text();
    if(m_PhysicalBoundaryConditionsMap.contains(name)) {
      m_PBC_current = m_PhysicalBoundaryConditionsMap[name];
      m_PBC_current.setIndex(0);
    }
    else EG_BUG;
    loadPhysicalValues();
  }
}

//==========================

void GuiEditBoundaryConditions::setupSolvers()
{
  m_multipagewidget_Solver = new MultiPageWidget(ui.tab_Solver);
  m_multipagewidget_Solver->setObjectName(QString::fromUtf8("m_multipagewidget_Solver"));
  ui.verticalLayout_Solver->addWidget(m_multipagewidget_Solver);

  QFileInfo fileinfo;
  fileinfo.setFile(":/resources/solvers/solvers.txt");
  QFile file(fileinfo.filePath());
  if (!file.exists()) {
    qDebug() << "ERROR: " << fileinfo.filePath() << " not found.";
    EG_BUG;
  }
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    qDebug() << "ERROR:  Failed to open file " << fileinfo.filePath();
    EG_BUG;
  }
  QTextStream text_stream(&file);
  QString intext = text_stream.readAll();
  file.close();

  ///\todo Create a special parser method for this so that it can be reused by other classes
  int idx = 0;
  QStringList page_list = intext.split("=");
  foreach(QString page, page_list) {
    QString title;
    QString section;
    QString binary;
    QVector <QString> files;
    QStringList variable_list = page.split(";");
    foreach(QString variable, variable_list) {
      QStringList name_value = variable.split(":");
      if (name_value[0].trimmed() == "title") title = name_value[1].trimmed();
      if (name_value[0].trimmed() == "section") section = name_value[1].trimmed();
      if (name_value[0].trimmed() == "binary") binary = name_value[1].trimmed();
      if (name_value[0].trimmed() == "files") {
        QStringList file_list = name_value[1].split(",");
        foreach(QString file, file_list) {
          files.push_back(":/resources/solvers/" + section + "/" + file.trimmed());
        }
      }
    }

    m_SolverBinary.push_back(binary);
    MultiPageWidgetPage* page = new MultiPageWidgetPage(files, section, m_multipagewidget_Solver);
    m_page_vector.push_back(page);
    m_multipagewidget_Solver->addPage((QWidget*)page);
    m_multipagewidget_Solver->setPageTitle(title, idx);

    idx++;
  }

  m_multipagewidget_Solver->setCurrentIndex(GuiMainWindow::pointer()->getXmlSection("solver/general/solver_type").toInt());
}

void GuiEditBoundaryConditions::saveSolverParameters()
{
  //Save solver parameters
  for (int i = 0; i < m_page_vector.size(); i++) {
    m_page_vector[i]->saveEgc();
  }
  QString solver_type;
  int solver_type_index = m_multipagewidget_Solver->currentIndex();
  solver_type.setNum(solver_type_index);
  GuiMainWindow::pointer()->setXmlSection("solver/general/solver_type", solver_type);
  GuiMainWindow::pointer()->setXmlSection("solver/general/solver_binary", m_SolverBinary[solver_type_index]);
}

//==========================

void GuiEditBoundaryConditions::addProcess()
{
  int row = ui.tableWidget_Processes->currentRow();
  QString host, weight;
  if(row<0) {
    row = 0;
    host = "host";
    weight = "1";
  }
  else {
    host = ui.tableWidget_Processes->item(row, 0)->text();
    weight = ui.tableWidget_Processes->item(row, 1)->text();
  }
  qDebug()<<"row="<<row;
  ui.tableWidget_Processes->insertRow(row);
  ui.tableWidget_Processes->setItem(row, 0, new QTableWidgetItem());
  ui.tableWidget_Processes->setItem(row, 1, new QTableWidgetItem());
  ui.tableWidget_Processes->item(row, 0)->setText(host);
  ui.tableWidget_Processes->item(row, 1)->setText(weight);
  ui.tableWidget_Processes->resizeColumnsToContents();
}

void GuiEditBoundaryConditions::deleteProcess()
{
  ui.tableWidget_Processes->removeRow(ui.tableWidget_Processes->currentRow());
}

void GuiEditBoundaryConditions::importHostFile()
{

}

void GuiEditBoundaryConditions::exportHostFile()
{

}

void GuiEditBoundaryConditions::loadMpiParameters()
{
  QString hostfile_txt = GuiMainWindow::pointer()->getXmlSection( "solver/general/host_weight_list" );
  stringToTable(hostfile_txt);
}

void GuiEditBoundaryConditions::saveMpiParameters()
{
  QString hostfile_txt = tableToString();
  GuiMainWindow::pointer()->setXmlSection( "solver/general/host_weight_list", hostfile_txt );
}

QString GuiEditBoundaryConditions::tableToString()
{
  QString hostfile_txt;
  for(int row = 0; row < ui.tableWidget_Processes->rowCount(); row++) {
    QString host, weight;
    host = ui.tableWidget_Processes->item(row, 0)->text();
    weight = ui.tableWidget_Processes->item(row, 1)->text();
    hostfile_txt += host + ":" + weight;
    if(row!=ui.tableWidget_Processes->rowCount()-1) hostfile_txt += ",";
  }
  return(hostfile_txt);
}

void GuiEditBoundaryConditions::stringToTable(QString hostfile_txt)
{
  QVector <QString> host;
  QVector <QString> weight;
  
  QStringList host_weight_list = hostfile_txt.split(",");
  foreach(QString host_weight, host_weight_list) {
    if(!host_weight.isEmpty()){
      QStringList values = host_weight.split(":");
//       qWarning()<<"values="<<values;
      host.push_back(values[0].trimmed());
      weight.push_back(values[1].trimmed());
    }
  }
  
  while (ui.tableWidget_Processes->rowCount()) {
    ui.tableWidget_Processes->removeRow(0);
  }
  
  for(int i = 0; i < host.size(); i++) {
    int row = ui.tableWidget_Processes->rowCount();
    ui.tableWidget_Processes->insertRow(row);
    ui.tableWidget_Processes->setItem(row, 0, new QTableWidgetItem());
    ui.tableWidget_Processes->setItem(row, 1, new QTableWidgetItem());
    ui.tableWidget_Processes->item(row, 0)->setText(host[i]);
    ui.tableWidget_Processes->item(row, 1)->setText(weight[i]);
    ui.tableWidget_Processes->resizeColumnsToContents();
  }
}
