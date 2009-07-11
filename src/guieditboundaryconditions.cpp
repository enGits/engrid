//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#include "physicalboundaryconditions.h"
#include "multipagewidget.h"
#include "multipagewidgetpage.h"

#include <QVBoxLayout>
#include <QFileInfo>

GuiEditBoundaryConditions::GuiEditBoundaryConditions()
{
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
  resetOrientation(grid);
  while (ui.T->rowCount()) ui.T->removeRow(0);
  foreach (int i, boundary_codes) {
    BoundaryCondition bc = (*m_BcMap)[i];
    ui.T->insertRow(ui.T->rowCount());
    int r = ui.T->rowCount()-1;
    ui.T->setItem(r,0,new QTableWidgetItem());
    ui.T->item(r,0)->setFlags(ui.T->item(r,0)->flags() & (~Qt::ItemIsSelectable));
    ui.T->item(r,0)->setFlags(ui.T->item(r,0)->flags() & (~Qt::ItemIsEditable));
    ui.T->setItem(r,1,new QTableWidgetItem());
    ui.T->setItem(r,2,new QTableWidgetItem());
    QString idx;
    idx.setNum(i);
    ui.T->item(r,0)->setText(idx);
    QString name = bc.getName();
    if (name == "unknown") name = QString("BC") + idx;
    ui.T->item(r,1)->setText(name);
    ui.T->item(r,2)->setText(bc.getType());
  }
  
  updateVol();
  updatePhysicalBoundaryConditions();
  
  m_PreviousSelected = 0;
  ui.listWidget_BoundaryType->setCurrentRow(m_PreviousSelected);
  if(ui.listWidget_BoundaryType->count()>0) loadPhysicalValues(ui.listWidget_BoundaryType->currentItem()->text());
  
  connect(ui.pushButtonAdd, SIGNAL(clicked()), this, SLOT(addVol()));
  connect(ui.pushButtonDelete, SIGNAL(clicked()), this, SLOT(delVol()));
  connect(ui.pushButton_AddBoundaryType, SIGNAL(clicked()), this, SLOT(addBoundaryType()));
  connect(ui.pushButton_DeleteBoundaryType, SIGNAL(clicked()), this, SLOT(deleteBoundaryType()));
  connect(ui.listWidget_BoundaryType, SIGNAL(itemSelectionChanged()), this, SLOT(changePhysicalValues()));
  
  ui.T->setItemDelegate(delegate);
}

void GuiEditBoundaryConditions::setupSolvers()
{
  m_multipagewidget_Solver = new MultiPageWidget(ui.tab_Solver);
  m_multipagewidget_Solver->setObjectName(QString::fromUtf8("m_multipagewidget_Solver"));
  ui.verticalLayout_Solver->addWidget(m_multipagewidget_Solver);
  
  QFileInfo fileinfo;
  fileinfo.setFile( ":/resources/solvers/solvers.txt" );
  QFile file( fileinfo.filePath() );
  if ( !file.exists() ) {
    qDebug() << "ERROR: " << fileinfo.filePath() << " not found.";
    EG_BUG;
  }
  if ( !file.open( QIODevice::ReadOnly | QIODevice::Text ) ) {
    qDebug() << "ERROR:  Failed to open file " << fileinfo.filePath();
    EG_BUG;
  }
  QTextStream text_stream( &file );
  QString intext = text_stream.readAll();
  file.close();
  
  int idx = 0;
  QStringList page_list = intext.split("=");
  foreach(QString page, page_list) {
    QString title;
    QString section;
    QVector <QString> files;
    QStringList variable_list = page.split(";");
    foreach(QString variable, variable_list) {
      QStringList name_value = variable.split(":");
      if(name_value[0].trimmed()=="title") title = name_value[1].trimmed();
      if(name_value[0].trimmed()=="section") section = name_value[1].trimmed();
      if(name_value[0].trimmed()=="files") {
        QStringList file_list = name_value[1].split(",");
        foreach(QString file, file_list) {
          files.push_back(":/resources/solvers/" + section + "/" + file.trimmed());
        }
      }
    }
    qDebug()<<"title="<<title;
    qDebug()<<"section="<<section;
    qDebug()<<"files="<<files;
    MultiPageWidgetPage* page = new MultiPageWidgetPage(files, section, m_multipagewidget_Solver);
    m_page_vector.push_back(page);
    m_multipagewidget_Solver->addPage( (QWidget*)page );
    m_multipagewidget_Solver->setPageTitle(title, idx);
    idx++;
  }
  
  m_multipagewidget_Solver->setCurrentIndex(GuiMainWindow::pointer()->getXmlSection("solver/general/solver_type").toInt());
}

void GuiEditBoundaryConditions::saveSolverParameters()
{
  //Save solver parameters
  for(int i = 0; i < m_page_vector.size(); i++) {
    m_page_vector[i]->saveEgc();
  }
  QString solver_type;
  solver_type.setNum(m_multipagewidget_Solver->currentIndex());
  GuiMainWindow::pointer()->setXmlSection("solver/general/solver_type", solver_type);
}

void GuiEditBoundaryConditions::loadPhysicalValues(QString name)
{
  if(m_PhysicalBoundaryConditionsMap.contains(name)) {
    PhysicalBoundaryConditions PBC = m_PhysicalBoundaryConditionsMap[name];
    QString str;
    str.setNum(PBC.m_Pressure); ui.lineEditPressure->setText(str);
    str.setNum(PBC.m_Temperature); ui.lineEditTemperature->setText(str);
    str.setNum(PBC.m_Velocity); ui.lineEditVelocity->setText(str);
  }
}

void GuiEditBoundaryConditions::savePhysicalValues(QString name, int index)
{
  if(m_PhysicalBoundaryConditionsMap.contains(name)) {
    PhysicalBoundaryConditions PBC(name, index);
    PBC.m_Pressure = ui.lineEditPressure->text().toDouble();
    PBC.m_Temperature = ui.lineEditTemperature->text().toDouble();
    PBC.m_Velocity = ui.lineEditVelocity->text().toDouble();
    m_PhysicalBoundaryConditionsMap[PBC.getName()] = PBC;
  }
}

void GuiEditBoundaryConditions::changePhysicalValues()
{
  if(ui.listWidget_BoundaryType->count()>0) {
    int index = ui.listWidget_BoundaryType->currentRow();
    QString name = ui.listWidget_BoundaryType->currentItem()->text();
    
    cout<<"switching to index = "<< index << " name = " << qPrintable(name) << endl;
    
    savePhysicalValues(m_PreviousSelectedName, m_PreviousSelectedIndex);
    loadPhysicalValues(name);
    
    m_PreviousSelectedName = name;
    m_PreviousSelectedIndex = index;
  }
}

void GuiEditBoundaryConditions::addBoundaryType()
{
  cout<<"Adding BT"<<endl;
  QString name = ui.lineEdit_BoundaryTypes->text();
  if (!name.isEmpty() && !m_PhysicalBoundaryConditionsMap.contains(name)) {
    PhysicalBoundaryConditions PBC(name, ui.listWidget_BoundaryType->count());
    m_PhysicalBoundaryConditionsMap[PBC.getName()] = PBC;
    ui.listWidget_BoundaryType->addItem(PBC.getName());
  }
}

void GuiEditBoundaryConditions::deleteBoundaryType()
{
  cout<<"Deleting BT"<<endl;
  int row = ui.listWidget_BoundaryType->currentRow();
  cout<<"row="<<row<<endl;
  if(ui.listWidget_BoundaryType->count()<=1) {
    ui.listWidget_BoundaryType->clear();
    m_PhysicalBoundaryConditionsMap.clear();
  }
  else if( 0 <= row && row < ui.listWidget_BoundaryType->count() ) {
    QListWidgetItem* list_widget_item = ui.listWidget_BoundaryType->takeItem(row);
    m_PhysicalBoundaryConditionsMap.remove(list_widget_item->text());
    delete list_widget_item;
  }
  else {
    cout<<"Nothing to delete."<<endl;
  }
}

void GuiEditBoundaryConditions::updatePhysicalBoundaryConditions()
{
  ui.listWidget_BoundaryType->clear();
  QList<PhysicalBoundaryConditions> physical_boundary_conditions = GuiMainWindow::pointer()->getAllPhysicalBoundaryConditions();
  foreach (PhysicalBoundaryConditions PBC, physical_boundary_conditions) {
    ui.listWidget_BoundaryType->addItem(PBC.getName());
    m_PhysicalBoundaryConditionsMap[PBC.getName()] = PBC;
  }
}

void GuiEditBoundaryConditions::updateVol()
{
  while (ui.T->columnCount() > 3) {
    ui.T->removeColumn(3);
  }
  QList<VolumeDefinition> vols = GuiMainWindow::pointer()->getAllVols();
  foreach (VolumeDefinition V, vols) {
    int c = ui.T->columnCount();
    ui.T->insertColumn(c);
    ui.T->setHorizontalHeaderItem(c, new QTableWidgetItem(V.getName()));
    for (int i = 0; i < ui.T->rowCount(); ++i) {
      int bc = ui.T->item(i,0)->text().toInt();
      if      (V.getSign(bc) == 1)  ui.T->setItem(i, c, new QTableWidgetItem("green"));
      else if (V.getSign(bc) == -1) ui.T->setItem(i, c, new QTableWidgetItem("yellow"));
      else                          ui.T->setItem(i, c, new QTableWidgetItem(" "));
    }
  }
}

void GuiEditBoundaryConditions::addVol()
{
  QString name = ui.lineEditVolume->text();
  if (!name.isEmpty()) {
    VolumeDefinition NV(name, ui.T->columnCount()-2);
    QList<VolumeDefinition> vols;
    QList<VolumeDefinition> new_vols;
    vols = GuiMainWindow::pointer()->getAllVols();
    foreach (VolumeDefinition V, vols) {
      if (NV.getName() != V.getName()) {
        new_vols.push_back(V);
      }
    }
    new_vols.push_back(NV);
    GuiMainWindow::pointer()->setAllVols(new_vols);
    updateVol();
  }
}

///@@@ TODO: Fix bug of reloading all volumes when deleting volumes and then adding new ones
void GuiEditBoundaryConditions::delVol()
{
  int c = ui.T->currentColumn();
  if (c > 2) {
    ui.T->removeColumn(c);
  }
}

void GuiEditBoundaryConditions::operate()
{
  // BoundaryCondition and VolumeDefinition
  QVector<VolumeDefinition> vols(ui.T->columnCount());
  for (int j = 3; j < ui.T->columnCount(); ++j) {
    QString vol_name = ui.T->horizontalHeaderItem(j)->text();
    VolumeDefinition V(vol_name, j-2);
    vols[j] = V;
  }
  for (int i = 0; i < ui.T->rowCount(); ++i) {
    int bc = ui.T->item(i,0)->text().toInt();
    BoundaryCondition BC(ui.T->item(i,1)->text(),ui.T->item(i,2)->text());
    (*m_BcMap)[bc] = BC;
    for (int j = 3; j < ui.T->columnCount(); ++j) {
      QString vol_name = ui.T->horizontalHeaderItem(j)->text();
      VolumeDefinition V = vols[j];
      if      (ui.T->item(i,j)->text() == "green")  V.addBC(bc,  1);
      else if (ui.T->item(i,j)->text() == "yellow") V.addBC(bc, -1);
      else                                          V.addBC(bc,  0);
      vols[j] = V;
    }
  }
  QList<VolumeDefinition> vol_list;
  for (int j = 3; j < ui.T->columnCount(); ++j) {
    vol_list.append(vols[j]);
  }
  GuiMainWindow::pointer()->setAllVols(vol_list);
  
  // PhysicalBoundaryConditions
  GuiMainWindow::pointer()->setAllPhysicalBoundaryConditions(m_PhysicalBoundaryConditionsMap);
  
  saveSolverParameters();
  saveMpiParameters();
}

void GuiEditBoundaryConditions::saveMpiParameters()
{
  GuiMainWindow::pointer()->setXmlSection("solver/general/hostfile", ui.plainTextEdit_HostFile->toPlainText());
  QString str_num_processes;
  str_num_processes.setNum(ui.spinBox_NumProcesses->value());
  GuiMainWindow::pointer()->setXmlSection("solver/general/num_processes", str_num_processes);
}

void GuiEditBoundaryConditions::loadMpiParameters()
{
  ui.plainTextEdit_HostFile->setPlainText(GuiMainWindow::pointer()->getXmlSection("solver/general/hostfile"));
  ui.spinBox_NumProcesses->setValue(GuiMainWindow::pointer()->getXmlSection("solver/general/num_processes").toInt());
}
