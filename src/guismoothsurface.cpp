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
#include "guismoothsurface.h"

#include "swaptriangles.h"
#include "surfacemesher.h"
#include "vertexdelegate.h"
#include "settingssheet.h"
#include "laplacesmoother.h"
#include "guimainwindow.h"
#include "containertricks.h"
#include "updatedesiredmeshdensity.h"

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

#include <stdio.h>

///////////////////////////////////////////
/* Here is how we we get QTextStreams that look like iostreams */
QTextStream Qcin(stdin, QIODevice::ReadOnly);
QTextStream Qcout(stdout, QIODevice::WriteOnly);
QTextStream Qcerr(stderr, QIODevice::WriteOnly);
///////////////////////////////////////////

//////////////////////////////////////////////

int cout_vtkWindowedSincPolyDataFilter(vtkWindowedSincPolyDataFilter* smooth)
{
  cout<<"NumberOfIterations="<<smooth->GetNumberOfIterations()<<endl;
  cout<<"PassBand="<<smooth->GetPassBand()<<endl;
  cout<<"FeatureEdgeSmoothing="<<smooth->GetFeatureEdgeSmoothing()<<endl;
  cout<<"FeatureAngle="<<smooth->GetFeatureAngle()<<endl;
  cout<<"EdgeAngle="<<smooth->GetEdgeAngle()<<endl;
  cout<<"BoundarySmoothing="<<smooth->GetBoundarySmoothing()<<endl;
  cout<<"GenerateErrorScalars="<<smooth->GetGenerateErrorScalars()<<endl;
  cout<<"GenerateErrorVectors="<<smooth->GetGenerateErrorVectors()<<endl;
  return(0);
}

int cout_vtkSmoothPolyDataFilter(vtkSmoothPolyDataFilter* smooth)
{
  cout<<"GetConvergence="<<smooth->GetConvergence()<<endl;
  cout<<"GetNumberOfIterations="<<smooth->GetNumberOfIterations()<<endl;
  cout<<"GetRelaxationFactor="<<smooth->GetRelaxationFactor()<<endl;
  cout<<"getFeatureEdgeSmoothing="<<smooth->GetFeatureEdgeSmoothing()<<endl;
  cout<<"GetFeatureAngle="<<smooth->GetFeatureAngle()<<endl;
  cout<<"GetEdgeAngle="<<smooth->GetEdgeAngle()<<endl;
  cout<<"GetBoundarySmoothing="<<smooth->GetBoundarySmoothing()<<endl;
  cout<<"GetGenerateErrorScalars="<<smooth->GetGenerateErrorScalars()<<endl;
  cout<<"GetGenerateErrorVectors="<<smooth->GetGenerateErrorVectors()<<endl;
  return(0);
}

GuiSmoothSurface::GuiSmoothSurface()
{
  setQuickSave(true);

  m_tableWidget = new SettingsSheet();
  ui.verticalLayout_SettingsSheet->addWidget(m_tableWidget);
  
  populateBoundaryCodes(ui.listWidget);
  populateBoundaryCodes(ui.listWidget_Source);
  
  ui.SmoothMethod->addItem("Method 0: vtkSmoothPolyDataFilter smoothing");
  ui.SmoothMethod->addItem("Method 1: vtkWindowedSincPolyDataFilter smoothing");
  ui.SmoothMethod->addItem("Method 2: Laplacian smoothing");
  ui.SmoothMethod->addItem("Method 3: swap triangles");
  ui.SmoothMethod->addItem("Method 4: Refine mesh");
  ui.SmoothMethod->addItem("Method 5: Update node information");
  ui.SmoothMethod->addItem("Method 6: Decimate (Delete all possible points)");
  ui.SmoothMethod->addItem("Method 7: Delete selected points");
  ui.SmoothMethod->addItem("Method 8: Projection test");
  ui.SmoothMethod->addItem("Method 9: Save selected boundary codes");
  
  vtkSmoothPolyDataFilter* smooth=vtkSmoothPolyDataFilter::New();
  vtkWindowedSincPolyDataFilter* smooth2=vtkWindowedSincPolyDataFilter::New();
  
  ui.doubleSpinBox_Convergence->setValue(smooth->GetConvergence());
  ui.spinBox_NumberOfIterations->setValue(smooth->GetNumberOfIterations());
  QString tmp;
  ui.lineEdit_RelaxationFactor->setText(tmp.setNum(smooth->GetRelaxationFactor()));
  ui.doubleSpinBox_PassBand->setValue(smooth2->GetPassBand());
  ui.checkBox_FeatureEdgeSmoothing->setCheckState(Qt::Checked);
  ui.doubleSpinBox_FeatureAngle->setValue(smooth->GetFeatureAngle());
  ui.doubleSpinBox_EdgeAngle->setValue(smooth->GetEdgeAngle());
  ui.checkBox_BoundarySmoothing->setCheckState(int2CheckState(Qt::Checked));
  ui.checkBox_GenerateErrorScalars->setCheckState(int2CheckState(smooth->GetGenerateErrorScalars()));
  ui.checkBox_GenerateErrorVectors->setCheckState(int2CheckState(smooth->GetGenerateErrorVectors()));
  
  smooth->Delete();
  smooth2->Delete();
  
  //Load settings
  readSettings();
  
  int Nrow,Ncol;
  Nrow=m_tableWidget->rowCount();
  Ncol=m_tableWidget->columnCount();
  
  QList<QString> list;
  list
    <<"VTK_SIMPLE_VERTEX"
    <<"VTK_FIXED_VERTEX"
    <<"VTK_FEATURE_EDGE_VERTEX"
    <<"VTK_BOUNDARY_EDGE_VERTEX"
    <<"any";
  
  QList<QString> list2;
  list2
    << "yes"
    <<"no"
    <<"any";
  
  Nbc=ui.listWidget-> count ();
  m_tableWidget->setColumnCount(Nbc+3);
  VertexDelegate* item_delegate = new VertexDelegate(Nbc, list);
  m_tableWidget->setItemDelegate(item_delegate);
  
  QStringList L;
  for(int i=0;i<Nbc;i++)
  {
    L<<ui.listWidget->item(i)->text();
  }
  L<<"Vertex Type";
  L<<"Nodelist";
  L<<"Mesh Density";
  m_tableWidget->setHorizontalHeaderLabels(L);
  m_tableWidget->resizeColumnsToContents();
  
  
  current_filename= GuiMainWindow::pointer()->getFilename();
  Qcout<<"current_filename="<<current_filename<<endl;
  Qcout<<"Loading settings from "+current_filename+".sp..."<<endl;
  
  if (!m_tableWidget->readFile(current_filename+".sp",0)) {
    cout<<"Loading settingssheet failed"<<endl;
  }
  
  connect(ui.pushButton_AddSet, SIGNAL(clicked()), this, SLOT(AddSet()));
  connect(ui.pushButton_RemoveSet, SIGNAL(clicked()), this, SLOT(RemoveSet()));
  connect(ui.pushButton_TestSet, SIGNAL(clicked()), this, SLOT(TestSet()));
  connect(ui.pushButton_Load, SIGNAL(clicked()), this, SLOT(Load()));
  connect(ui.pushButton_Save, SIGNAL(clicked()), this, SLOT(Save()));
  connect(ui.pushButton_SelectAll_BC, SIGNAL(clicked()), this, SLOT(SelectAll_BC()));
  connect(ui.pushButton_ClearAll_BC, SIGNAL(clicked()), this, SLOT(ClearAll_BC()));
  connect(ui.pushButton_SelectAll_Source, SIGNAL(clicked()), this, SLOT(SelectAll_Source()));
  connect(ui.pushButton_ClearAll_Source, SIGNAL(clicked()), this, SLOT(ClearAll_Source()));
}

///////////////////////////////////////////

int GuiSmoothSurface::readSettings()
{
  QSettings local_qset("enGits","enGrid_smoothsurface");
//   current_settingssheet_name=local_qset.value("Filename", "").toString();
  ui.SmoothMethod->setCurrentIndex(local_qset.value("Method", 0).toInt());
  ui.spinBox_NumberOfSmoothIterations->setValue(local_qset.value("NumberOfSmoothIterations", 20).toInt());
  ui.spinBox_maxiter_density->setValue(local_qset.value("maxiter_density", 1000).toInt());
  ui.spinBox_DebugLevel->setValue(local_qset.value("DebugLevel", 0).toInt());
  ui.spinBox_NumberOfIterations->setValue(local_qset.value("NumberOfIterations", 1).toInt());
  
  ui.doubleSpinBox_Convergence_meshdensity->setValue(local_qset.value("Convergence_meshdensity", 0.000001).toDouble());
  ui.checkBox_insert_FP->setCheckState(int2CheckState(local_qset.value("insert_FP", 2).toInt()));
  ui.checkBox_insert_EP->setCheckState(int2CheckState(local_qset.value("insert_EP", 2).toInt()));
  ui.checkBox_remove_FP->setCheckState(int2CheckState(local_qset.value("remove_FP", 2).toInt()));
  ui.checkBox_remove_EP->setCheckState(int2CheckState(local_qset.value("remove_EP", 2).toInt()));
  
  ui.checkBox_GenerateErrorScalars->setCheckState(int2CheckState(local_qset.value("GenerateErrorScalars", 2).toInt()));
  ui.checkBox_GenerateErrorVectors->setCheckState(int2CheckState(local_qset.value("GenerateErrorVectors", 2).toInt()));
  
  ui.checkBox_Swap->setCheckState(int2CheckState(local_qset.value("DoSwap", 2).toInt()));
  ui.checkBox_LaplaceSmoothing->setCheckState(int2CheckState(local_qset.value("DoLaplaceSmoothing", 2).toInt()));
  
  if(local_qset.value("DensityUnit_is_length", false).toBool()){
    ui.radioButton_length->toggle();
  }
  else{
    ui.radioButton_density->toggle();
  }
  
  int size;
  size = local_qset.beginReadArray("list_BC");
  if(ui.listWidget->count()==size)
  {
    for (int i = 0; i < size; ++i) {
      local_qset.setArrayIndex(i);
      Qt::CheckState x=int2CheckState(local_qset.value("state").toInt());
      ui.listWidget->item(i)->setCheckState(x);
    }
    local_qset.endArray();
  }
  
  size = local_qset.beginReadArray("list_BC_Source");
  if(ui.listWidget_Source->count()==size)
  {
    for (int i = 0; i < size; ++i) {
      local_qset.setArrayIndex(i);
      Qt::CheckState x=int2CheckState(local_qset.value("state").toInt());
      ui.listWidget_Source->item(i)->setCheckState(x);
    }
    local_qset.endArray();
  }
  return(0);
}

int GuiSmoothSurface::writeSettings()
{
  QSettings local_qset("enGits","enGrid_smoothsurface");
//   local_qset.setValue("Filename", current_settingssheet_name);
  local_qset.setValue("Method", ui.SmoothMethod->currentIndex());
  local_qset.setValue("NumberOfSmoothIterations", ui.spinBox_NumberOfSmoothIterations->value());
  local_qset.setValue("NumberOfIterations", ui.spinBox_NumberOfIterations->value());
  local_qset.setValue("maxiter_density", ui.spinBox_maxiter_density->value());
  local_qset.setValue("DebugLevel", ui.spinBox_DebugLevel->value());
  local_qset.setValue("Convergence_meshdensity", ui.doubleSpinBox_Convergence_meshdensity->value());
  
  local_qset.setValue("insert_FP", ui.checkBox_insert_FP->checkState());
  local_qset.setValue("insert_EP", ui.checkBox_insert_EP->checkState());
  local_qset.setValue("remove_FP", ui.checkBox_remove_FP->checkState());
  local_qset.setValue("remove_EP", ui.checkBox_remove_EP->checkState());
  
  local_qset.setValue("GenerateErrorScalars", ui.checkBox_GenerateErrorScalars->checkState());
  local_qset.setValue("GenerateErrorVectors", ui.checkBox_GenerateErrorVectors->checkState());
  
  local_qset.setValue("DoSwap", ui.checkBox_Swap->checkState());
  local_qset.setValue("DoLaplaceSmoothing", ui.checkBox_LaplaceSmoothing->checkState());
  local_qset.setValue("DensityUnit_is_length",ui.radioButton_length->isChecked());
  
  QList<Qt::CheckState> list;
  
  for (int i = 0; i < ui.listWidget->count(); ++i) {
    list << ui.listWidget->item(i)->checkState();
  };
  local_qset.beginWriteArray("list_BC");
  for (int i = 0; i < list.size(); ++i) {
    local_qset.setArrayIndex(i);
    local_qset.setValue("state", list.at(i));
  }
  local_qset.endArray();
  
  list.clear();
  for (int i = 0; i < ui.listWidget_Source->count(); ++i) {
    list << ui.listWidget_Source->item(i)->checkState();
  };
  local_qset.beginWriteArray("list_BC_Source");
  for (int i = 0; i < list.size(); ++i) {
    local_qset.setArrayIndex(i);
    local_qset.setValue("state", list.at(i));
  }
  local_qset.endArray();
  return(0);
}

///////////////////////////////////////////

void GuiSmoothSurface::before()
{

}

void GuiSmoothSurface::Load()
{
  QString current_settingssheet_name = QFileDialog::getOpenFileName(this,tr("Open SettingsSheet"), ".",tr("SettingsSheet files (*.sp)"));
  if (!current_settingssheet_name.isEmpty() && !m_tableWidget->readFile(current_settingssheet_name)) {
    cout<<"Loading failed"<<endl;
  }
}
void GuiSmoothSurface::Save()
{
  QString current_settingssheet_name = QFileDialog::getSaveFileName(this,tr("Save SettingsSheet as..."), ".",tr("SettingsSheet files (*.sp)"));
  if (!current_settingssheet_name.isEmpty() && !m_tableWidget->writeFile(current_settingssheet_name)) {
    cout<<"Saving failed"<<endl;
  }
}

void GuiSmoothSurface::SelectAll_BC()
{
  for (int i = 0; i < ui.listWidget->count(); ++i) {
    ui.listWidget->item(i)->setCheckState(Qt::Checked);
  };
}
void GuiSmoothSurface::ClearAll_BC()
{
  for (int i = 0; i < ui.listWidget->count(); ++i) {
    ui.listWidget->item(i)->setCheckState(Qt::Unchecked);
  };
}
void GuiSmoothSurface::SelectAll_Source()
{
  for (int i = 0; i < ui.listWidget_Source->count(); ++i) {
    ui.listWidget_Source->item(i)->setCheckState(Qt::Checked);
  };
}
void GuiSmoothSurface::ClearAll_Source()
{
  for (int i = 0; i < ui.listWidget_Source->count(); ++i) {
    ui.listWidget_Source->item(i)->setCheckState(Qt::Unchecked);
  };
}

void GuiSmoothSurface::TestSet()
{
  cout<<"Testing set"<<endl;
  getSet();
}

//This is where we get the user defined mesh densities
QVector <VertexMeshDensity> GuiSmoothSurface::getSet()
{
  cout<<"Getting set"<<endl;
  QVector <VertexMeshDensity> VMDvector;
  
  cout<<"VMDvector:"<<VMDvector<<endl;
  
  int N_VMD=m_tableWidget->rowCount();
  VMDvector.resize(N_VMD);
  cout<<"VMDvector.size()="<<VMDvector.size()<<endl;
  for(int i=0;i<N_VMD;i++)
  {
    for(int j=0;j<Nbc;j++)
    {
      int bc = m_tableWidget->horizontalHeaderItem(j)->text().toInt();
      int state = CheckState2int( m_tableWidget->item(i,j)->checkState() );
      VMDvector[i].BCmap[bc]=state;
    }
    VMDvector[i].type=Str2VertexType(m_tableWidget->item(i,Nbc)->text());
    VMDvector[i].setNodes(m_tableWidget->item(i,Nbc+1)->text());
    if(ui.radioButton_density->isChecked()){
      VMDvector[i].density=m_tableWidget->item(i,Nbc+2)->text().toDouble();
    }
    else{
      cout<<"desired_density="<<1.0/(m_tableWidget->item(i,Nbc+2)->text().toDouble())<<endl;
      VMDvector[i].density=1.0/(m_tableWidget->item(i,Nbc+2)->text().toDouble());
    }
  }
  cout<<"VMDvector:"<<VMDvector<<endl;
  return(VMDvector);
}

void GuiSmoothSurface::AddSet()
{
  cout<<"Adding set"<<endl;
  int row=m_tableWidget->rowCount();
  m_tableWidget->insertRow(row);
  
  int Nbc=ui.listWidget->count();
  for(int i=0;i<Nbc;i++)
  {
    TriStateTableWidgetItem* newBC = new TriStateTableWidgetItem();
    newBC->setFlags(Qt::ItemIsTristate | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
    m_tableWidget->setItem(row, i, newBC);
  }
  QTableWidgetItem* item1 = new QTableWidgetItem("any");
  m_tableWidget->setItem(row, Nbc, item1);
  QTableWidgetItem* item2 = new QTableWidgetItem("");
  m_tableWidget->setItem(row, Nbc+1, item2);
  QTableWidgetItem* item3 = new QTableWidgetItem("-1");
  m_tableWidget->setItem(row, Nbc+2, item3);
  m_tableWidget->resizeColumnsToContents();
}

void GuiSmoothSurface::RemoveSet()
{
  cout<<"Removing set"<<endl;
  m_tableWidget->removeRow(m_tableWidget->currentRow());
  m_tableWidget->resizeColumnsToContents();
}

int GuiSmoothSurface::DisplayErrorScalars(vtkPolyDataAlgorithm* algo)
{

  cout<<"==============="<<endl;
  cout<<"ErrorScalars:"<<endl;
  int N1,N2;
  double dist;
  N1=algo->GetOutput()->GetPointData()->GetNumberOfArrays();
//   cout<<"nb of arrays="<<N1<<endl;
  algo->GetOutput();//vtkPolyData
//   cout<<algo->GetOutput()->GetPointData()<<endl;//vtkPointData*
  vtkFloatArray *newScalars = vtkFloatArray::New();
  newScalars=(vtkFloatArray *)algo->GetOutput()->GetPointData()->GetArray(1);
  N1=newScalars->GetNumberOfComponents();
  N2=newScalars->GetNumberOfTuples();
  cout<<"Number of components=N1="<<N1<<endl;
  cout<<"Number of tuples=N2="<<N2<<endl;
  for (int i=0; i<N2; i++)
  {
    dist=newScalars->GetComponent(i-1,1);//strange, but works. O.o
    cout<<"dist["<<i<<"]="<<dist<<endl;
  }
  
  cout<<"==============="<<endl;
  newScalars->Delete();
  return(0);
}

int GuiSmoothSurface::DisplayErrorVectors(vtkPolyDataAlgorithm* algo)
{
  cout<<"==============="<<endl;
  cout<<"ErrorVectors:"<<endl;
//   int N1,N2;
/*  vec3_t xx;
  algo->GetOutput()->GetPoint(0, xx.data());*/
  int N1=algo->GetOutput()->GetPointData()->GetVectors()->GetNumberOfComponents();
  int N2=algo->GetOutput()->GetPointData()->GetVectors()->GetNumberOfTuples();
  cout<<"Number of components=N1="<<N1<<endl;
  cout<<"Number of tuples=N2="<<N2<<endl;
  
  ///@@@  TODO: Fix this eventually
/*  void vtkFieldData::GetTuple  	(  	const vtkIdType   	 i,
                               	   	double *  	tuple	 
                               	) 			
    Copy the ith tuple value into a user provided tuple array. Make sure that you've allocated enough space for the copy.
      Deprecated:
      as of VTK 5.2. Using this method for FieldData having arrays that are not subclasses of vtkDataArray may yield unexpected results. */
  //Yes, indeed, very unexpected... And what should we do instead?
  
  for(vtkIdType i=0;i<N2;i++)
  {
    double tuple[4];
    algo->GetOutput()->GetPointData()->GetTuple(i,tuple);
    ///@@@  TODO: This works, but seems incorrect
    cout<<"tuple["<<tuple[0]<<"]=("<<tuple[1]<<","<<tuple[2]<<","<<tuple[3]<<")"<<endl;
  }
  cout<<"==============="<<endl;
  return(0);
}

void GuiSmoothSurface::operate()
{
  cout<<"Saving settings..."<<endl;
  //Save settings
  writeSettings();
  
  Qcout<<"current_filename="<<current_filename<<endl;
  Qcout<<"Saving settings as "+current_filename+".sp..."<<endl;
  if(!m_tableWidget->writeFile(current_filename+".sp")) cout<<"Saving settingssheet failed."<<endl;
  
  cout<<"METHOD "<<ui.SmoothMethod->currentIndex()<<endl;
  //can't use switch case because dynamic variables seem to be forbidden inside case statements
  //////////////////////////////////////////////////////////////////////////////////////////////
  if(ui.SmoothMethod->currentIndex()==0)//vtkSmoothPolyDataFilter smoothing
  {
    //preparations
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QVector<vtkIdType> cells;
    getSurfaceCells(bcs, cells, this->grid);
    EG_VTKSP(vtkPolyData, pdata);
    addToPolyData(cells, pdata, this->grid);
    EG_VTKSP(vtkSmoothPolyDataFilter, smooth);
    
    cout_vtkSmoothPolyDataFilter(smooth);
    
    //configure vtkSmoothPolyDataFilter
    smooth->SetInput(pdata);
    
    smooth->SetConvergence (ui.doubleSpinBox_Convergence->value());
    smooth->SetNumberOfIterations (ui.spinBox_NumberOfIterations->value());
    smooth->SetRelaxationFactor (ui.lineEdit_RelaxationFactor->text().toDouble());
    smooth->SetFeatureEdgeSmoothing (ui.checkBox_FeatureEdgeSmoothing->checkState());
    smooth->SetFeatureAngle (ui.doubleSpinBox_FeatureAngle->value());
    smooth->SetEdgeAngle (ui.doubleSpinBox_EdgeAngle->value());
    smooth->SetBoundarySmoothing (ui.checkBox_BoundarySmoothing->checkState());
    smooth->SetGenerateErrorScalars (ui.checkBox_GenerateErrorScalars->checkState());
    smooth->SetGenerateErrorVectors (ui.checkBox_GenerateErrorVectors->checkState());
    
    QSet<int> bcs_Source;
    getSelectedItems(ui.listWidget,bcs_Source);
    QVector<vtkIdType> cells_Source;
    getSurfaceCells(bcs_Source, cells_Source, this->grid);
    EG_VTKSP(vtkPolyData, pdata_Source);
    addToPolyData(cells_Source, pdata_Source, this->grid);
    smooth->SetSource (pdata_Source);
    
    cout_vtkSmoothPolyDataFilter(smooth);
    
    //smooth
    smooth->Update();
    
    if(ui.checkBox_GenerateErrorScalars->checkState()) DisplayErrorScalars(smooth);
    if(ui.checkBox_GenerateErrorVectors->checkState()) DisplayErrorVectors(smooth);
    
    //copy smoothed grid to main grid
    EG_VTKDCN(vtkLongArray_t, node_index, pdata, "node_index");
    for (vtkIdType i = 0; i < smooth->GetOutput()->GetNumberOfPoints(); ++i) {
      vec3_t x;
      smooth->GetOutput()->GetPoints()->GetPoint(i, x.data());
      vtkIdType nodeId = node_index->GetValue(i);
      this->grid->GetPoints()->SetPoint(nodeId, x.data());
    };
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==1)//vtkWindowedSincPolyDataFilter smoothing
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QVector<vtkIdType> cells;
    getSurfaceCells(bcs, cells, this->grid);
    EG_VTKSP(vtkPolyData, pdata);
    addToPolyData(cells, pdata, this->grid);
    EG_VTKSP(vtkWindowedSincPolyDataFilter, smooth);
    
    cout_vtkWindowedSincPolyDataFilter(smooth);
    
    smooth->SetInput(pdata);
    
    smooth->SetNumberOfIterations(ui.spinBox_NumberOfIterations->value());
    smooth->SetPassBand(ui.doubleSpinBox_PassBand->value());
    smooth->SetFeatureEdgeSmoothing(ui.checkBox_FeatureEdgeSmoothing->checkState());
    smooth->SetFeatureAngle(ui.doubleSpinBox_FeatureAngle->value());
    smooth->SetEdgeAngle(ui.doubleSpinBox_EdgeAngle->value());
    smooth->SetBoundarySmoothing(ui.checkBox_BoundarySmoothing->checkState());
    smooth->SetGenerateErrorScalars(ui.checkBox_GenerateErrorScalars->checkState());
    smooth->SetGenerateErrorVectors(ui.checkBox_GenerateErrorVectors->checkState());
    
    cout_vtkWindowedSincPolyDataFilter(smooth);
    
    smooth->Update();
    
    if(ui.checkBox_GenerateErrorScalars->checkState()) DisplayErrorScalars(smooth);
    if(ui.checkBox_GenerateErrorVectors->checkState()) DisplayErrorVectors(smooth);
    
    EG_VTKDCN(vtkLongArray_t, node_index, pdata, "node_index");
    for (vtkIdType i = 0; i < smooth->GetOutput()->GetNumberOfPoints(); ++i) {
      vec3_t x;
      smooth->GetOutput()->GetPoints()->GetPoint(i, x.data());
      vtkIdType nodeId = node_index->GetValue(i);
      this->grid->GetPoints()->SetPoint(nodeId, x.data());
    };
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==2)//Laplacian smoothing
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    
    LaplaceSmoother Lap;
    Lap.setGrid(this->grid);
    Lap.setBoundaryCodes(bcs);
    Lap.setSource(this->grid);
    Lap.setNumberOfIterations(ui.spinBox_NumberOfSmoothIterations->value());
    setDebugLevel(ui.spinBox_DebugLevel->value());
    Lap();
    Lap.delete_CellLocator_and_ProjectionSurface();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==3)//swap triangles
  {
    cout_grid(cout,this->grid);
    
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    
    QSet<int> bcs_complement=complementary_bcs(bcs,this->grid,cells);
    
    cout<<"bcs="<<bcs<<endl;
    cout<<"bcs_complement="<<bcs_complement<<endl;
    
    SwapTriangles swap;
    swap.setRespectBC(true);
    swap.setFeatureSwap(true);
    swap.setGrid(this->grid);
    swap.setBoundaryCodes(bcs_complement);
    swap();
    
    cout_grid(cout,this->grid);
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==4)// super smoothing
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    
    QVector <VertexMeshDensity> VMDvector=getSet();
    
    SurfaceMesher surfacemesher;
    
    surfacemesher.setGrid(this->grid);
    surfacemesher.setBoundaryCodes(bcs);
    
    surfacemesher.setVertexMeshDensityVector(VMDvector);
    surfacemesher.setConvergence (ui.doubleSpinBox_Convergence->value());
    surfacemesher.setNumberOfIterations (ui.spinBox_NumberOfIterations->value());
    surfacemesher.setRelaxationFactor (ui.lineEdit_RelaxationFactor->text().toDouble());
    surfacemesher.setFeatureEdgeSmoothing (ui.checkBox_FeatureEdgeSmoothing->checkState());
    surfacemesher.setFeatureAngle (ui.doubleSpinBox_FeatureAngle->value());
    surfacemesher.setEdgeAngle (ui.doubleSpinBox_EdgeAngle->value());
    surfacemesher.setBoundarySmoothing (ui.checkBox_BoundarySmoothing->checkState());
    
    surfacemesher.setConvergence_meshdensity(ui.doubleSpinBox_Convergence_meshdensity->value());
    
    surfacemesher.set_insert_FP(ui.checkBox_insert_FP->checkState());
    surfacemesher.set_insert_EP(ui.checkBox_insert_EP->checkState());
    surfacemesher.set_remove_FP(ui.checkBox_remove_FP->checkState());
    surfacemesher.set_remove_EP(ui.checkBox_remove_EP->checkState());
    surfacemesher.setDoSwap(ui.checkBox_Swap->checkState());
    surfacemesher.setDoLaplaceSmoothing(ui.checkBox_LaplaceSmoothing->checkState());
    
    surfacemesher.setN_SmoothIterations(ui.spinBox_NumberOfSmoothIterations->value());
    surfacemesher.setMaxiterDensity(ui.spinBox_maxiter_density->value());
    surfacemesher.setDebugLevel(ui.spinBox_DebugLevel->value());
    
    surfacemesher.setSource(this->grid);
    
    surfacemesher();
    surfacemesher.delete_CellLocator_and_ProjectionSurface();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==5)// Update current mesh density + node types + desired mesh density
  {
    setDebugLevel(ui.spinBox_DebugLevel->value());
    
    setConvergence(ui.doubleSpinBox_Convergence->value());
    setFeatureEdgeSmoothing(ui.checkBox_FeatureEdgeSmoothing->checkState());
    setFeatureAngle(ui.doubleSpinBox_FeatureAngle->value());
    setEdgeAngle(ui.doubleSpinBox_EdgeAngle->value());
    setBoundarySmoothing(ui.checkBox_BoundarySmoothing->checkState());
    
    UpdateCurrentMeshDensity();
    UpdateNodeType();
    
    QVector <VertexMeshDensity> VMDvector=getSet();
    
    UpdateDesiredMeshDensity update_desired_mesh_density;
    update_desired_mesh_density.setGrid(grid);
    update_desired_mesh_density.setConvergence_meshdensity(ui.doubleSpinBox_Convergence_meshdensity->value());
    update_desired_mesh_density.setMaxiterDensity(ui.spinBox_maxiter_density->value());
    update_desired_mesh_density.setVertexMeshDensityVector(VMDvector);
    update_desired_mesh_density();
    
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==6)// Delete all possible points
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    SurfaceMesher toto;
    
    setConvergence(ui.doubleSpinBox_Convergence->value());
    setFeatureEdgeSmoothing(ui.checkBox_FeatureEdgeSmoothing->checkState());
    setFeatureAngle(ui.doubleSpinBox_FeatureAngle->value());
    setEdgeAngle(ui.doubleSpinBox_EdgeAngle->value());
    setBoundarySmoothing(ui.checkBox_BoundarySmoothing->checkState());
    
    int N_newpoints;
    int N_newcells;
    
    bool Global_DelResult=true;
    while(Global_DelResult)
    {
      Global_DelResult=false;
      vtkIdType DeadNode=0;
      while(DeadNode<this->grid->GetNumberOfPoints())
      {
        bool Local_DelResult=true;
        while(Local_DelResult)
        {
          Local_DelResult = DeletePoint(DeadNode,N_newpoints,N_newcells);
          if(Local_DelResult) Global_DelResult=true;
        }
        DeadNode++;
      }
    }
    
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==7)// Delete selected points
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    SurfaceMesher toto;
    
    setDebugLevel(ui.spinBox_DebugLevel->value());
    
    setConvergence(ui.doubleSpinBox_Convergence->value());
    setFeatureEdgeSmoothing(ui.checkBox_FeatureEdgeSmoothing->checkState());
    setFeatureAngle(ui.doubleSpinBox_FeatureAngle->value());
    setEdgeAngle(ui.doubleSpinBox_EdgeAngle->value());
    setBoundarySmoothing(ui.checkBox_BoundarySmoothing->checkState());
    
    QVector <VertexMeshDensity> VMDvector=getSet();
    for(int i=0;i<VMDvector.size();i++)
    {
      cout<<"VMDvector["<<i<<"].nodeset="<<VMDvector[i].nodeset<<endl;
      int N_newpoints=0;
      int N_newcells=0;
      DeleteSetOfPoints(VMDvector[i].nodeset, N_newpoints, N_newcells);
      cout<<"N_newpoints="<<N_newpoints<<endl;
      cout<<"N_newcells="<<N_newcells<<endl;
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==8)// Projection test
  {
    //What we project on
    QSet<int> bcs_Source;
    getSelectedItems(ui.listWidget_Source, bcs_Source);
    
    //What we project
    QSet<int> bcs_Dest;
    getSelectedItems(ui.listWidget, bcs_Dest);
    
    EG_VTKSP(vtkUnstructuredGrid,grid_Source);
    getSurfaceCells(bcs_Source, cells, grid);
    getSubGrid(grid,cells,grid_Source);
    writeCells(grid,cells,GuiMainWindow::pointer()->getFilePath()+"Source.vtu");
    this->setSource(grid_Source);
    
    EG_VTKSP(vtkUnstructuredGrid,grid_Dest);
    makeCopy(grid,grid_Dest);
    getSurfaceCells(bcs_Dest, cells, grid_Dest);
    setCells(cells);
    
    foreach(vtkIdType id_node,nodes) {
      vec3_t M;
      grid_Dest->GetPoint(id_node,M.data());
      M = project(M);
      grid_Dest->GetPoints()->SetPoint(id_node,M.data());
    }
    
    makeCopy(grid_Dest,grid);
    
    this->delete_CellLocator_and_ProjectionSurface();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==9)// Save selected boundary codes
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    
    getSurfaceCells(bcs, cells, grid);
    writeCells(grid,cells,GuiMainWindow::pointer()->getFilePath()+"Selection.vtu");
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else
  {
    cout<<"UNKNOWN METHOD"<<endl;
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QVector<vtkIdType> cells;
    getSurfaceCells(bcs, cells, this->grid);
    setCells(cells);
    cout<<"cells="<<cells<<endl;
    cout<<"_cells="<<_cells<<endl;
    cout<<"nodes="<<nodes<<endl;
    cout<<"_nodes="<<_nodes<<endl;
    cout<<"n2c="<<n2c<<endl;
    cout<<"n2n="<<n2n<<endl;
    cout<<"c2c="<<c2c<<endl;
  }
};
