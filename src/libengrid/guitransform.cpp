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
#include "guitransform.h"
#include "seedsimpleprismaticlayer.h"
#include "gridsmoother.h"
#include "createvolumemesh.h"
#include "swaptriangles.h"
#include "deletetetras.h"
#include "deletecells.h"
#include <cmath>
#include "geometrytools.h"

using namespace GeometryTools;

/** Set default values */
void GuiTransform::before()
{
  cout<<"======================================"<<endl;
  cout<<"void GuiTransform::before()"<<endl;
  cout<<"======================================"<<endl;
  
  m_Ui.tabWidget->setCurrentIndex(0);
  
  m_Ui.lineEdit_Translation_X->setText("0");
  m_Ui.lineEdit_Translation_Y->setText("0");
  m_Ui.lineEdit_Translation_Z->setText("0");
  
  m_Ui.lineEdit_Rotation_Origin_X->setText("0");
  m_Ui.lineEdit_Rotation_Origin_Y->setText("0");
  m_Ui.lineEdit_Rotation_Origin_Z->setText("0");
  m_Ui.lineEdit_Rotation_Direction_X->setText("1");
  m_Ui.lineEdit_Rotation_Direction_Y->setText("0");
  m_Ui.lineEdit_Rotation_Direction_Z->setText("0");
  m_Ui.lineEdit_Rotation_Angle->setText("0");
  m_Ui.AngleInDegrees->setCheckState(Qt::Checked);
  
  m_Ui.lineEdit_Scaling_X->setText("1");
  m_Ui.lineEdit_Scaling_Y->setText("1");
  m_Ui.lineEdit_Scaling_Z->setText("1");
}

/** Apply transformations to the grid ( scaling, translation, rotation ) */
void GuiTransform::operate()
{
  l2g_t nodes = getPartNodes();

  cout<<"======================================"<<endl;
  cout<<"void GuiTransform::operate()"<<endl;
  cout<<"======================================"<<endl;
  
  //Translation
  vec3_t Translation_Vector( ( m_Ui.lineEdit_Translation_X->text()).toDouble(),
                             ( m_Ui.lineEdit_Translation_Y->text()).toDouble(),
                             ( m_Ui.lineEdit_Translation_Z->text()).toDouble());
  
  //Rotation
  vec3_t Rotation_Origin_Vector( ( m_Ui.lineEdit_Rotation_Origin_X->text()).toDouble(),
                                 ( m_Ui.lineEdit_Rotation_Origin_Y->text()).toDouble(),
                                 ( m_Ui.lineEdit_Rotation_Origin_Z->text()).toDouble());
  vec3_t Rotation_Direction_Vector( ( m_Ui.lineEdit_Rotation_Direction_X->text()).toDouble(),
                             ( m_Ui.lineEdit_Rotation_Direction_Y->text()).toDouble(),
                             ( m_Ui.lineEdit_Rotation_Direction_Z->text()).toDouble());
  double Rotation_Angle = ( m_Ui.lineEdit_Rotation_Angle->text()).toDouble();
  if(m_Ui.AngleInDegrees->checkState()) Rotation_Angle=deg2rad(Rotation_Angle);
  
  //Scaling
  vec3_t Scaling_Vector( ( m_Ui.lineEdit_Scaling_X->text()).toDouble(),
                             ( m_Ui.lineEdit_Scaling_Y->text()).toDouble(),
                             ( m_Ui.lineEdit_Scaling_Z->text()).toDouble());
  
  cout << "nodes.size()=" << nodes.size() << endl;
  cout << "Translation_Vector=" << Translation_Vector << endl;
  cout << "Rotation_Origin_Vector=" << Rotation_Origin_Vector << endl;
  cout << "Rotation_Direction_Vector=" << Rotation_Direction_Vector << endl;
  cout << "Rotation_Angle=" << Rotation_Angle << endl;
  cout << "AngleInDegrees=" << m_Ui.AngleInDegrees->checkState() << endl;
  cout << "Scaling_Vector=" << Scaling_Vector << endl;
  
  foreach(vtkIdType id_node, nodes)
  {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    
    x[0]+=Translation_Vector[0];
    x[1]+=Translation_Vector[1];
    x[2]+=Translation_Vector[2];
    
    Rotation_Direction_Vector.normalise();
    x = Rotation_Origin_Vector + GeometryTools::rotate(x-Rotation_Origin_Vector,Rotation_Direction_Vector,Rotation_Angle);
    
    x[0]*=Scaling_Vector[0];
    x[1]*=Scaling_Vector[1];
    x[2]*=Scaling_Vector[2];
    
    m_Grid->GetPoints()->SetPoint(id_node, x.data());
  }
  m_Grid->Modified();// to force a drawing update
};
