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

#include "guiedgelengthsourcebox.h"

GuiEdgeLengthSourceBox::GuiEdgeLengthSourceBox()
{
  m_X1 = vec3_t(0,0,0);
  m_X2 = vec3_t(1,0,0);
  m_Length = 1;
}

QString GuiEdgeLengthSourceBox::write()
{  
  QString txt = "box: " + name() + "; ";
  QString num;
  num.setNum(m_X1[0]); txt += num + "; ";
  num.setNum(m_X1[1]); txt += num + "; ";
  num.setNum(m_X1[2]); txt += num + "; ";
  num.setNum(m_X2[0]); txt += num + "; ";
  num.setNum(m_X2[1]); txt += num + "; ";
  num.setNum(m_X2[2]); txt += num + "; ";
  num.setNum(m_Length); txt += num + "; ";
  return txt;
}

bool GuiEdgeLengthSourceBox::read(QString txt)
{
  QStringList parts = txt.split(":");
  if (parts.size() != 2) {
    return false;
  }
  if (parts[0].trimmed() == "box") {
    QStringList words = parts[1].split(";");
    if (words.size() < 8) {
      return false;
    }
    m_Name = words[0].trimmed();
    m_X1 = vec3_t(words[1].toDouble(), words[2].toDouble(), words[3].toDouble());
    m_X2 = vec3_t(words[4].toDouble(), words[5].toDouble(), words[6].toDouble());
    m_Length = words[7].toDouble();
  } else {
    return false;
  }
  return true;
}

void GuiEdgeLengthSourceBox::setDlgFields()
{
  ui().lineEditName->setText(name());
  QString num;
  num.setNum(m_X1[0]); ui().lineEditX1->setText(num);
  num.setNum(m_X1[1]); ui().lineEditY1->setText(num);
  num.setNum(m_X1[2]); ui().lineEditZ1->setText(num);
  num.setNum(m_X2[0]); ui().lineEditX2->setText(num);
  num.setNum(m_X2[1]); ui().lineEditY2->setText(num);
  num.setNum(m_X2[2]); ui().lineEditZ2->setText(num);
  num.setNum(m_Length); ui().lineEditEL->setText(num);
}

void GuiEdgeLengthSourceBox::readDlgFields()
{
  setName(ui().lineEditName->text());
  m_X1[0]  = ui().lineEditX1->text().toDouble();
  m_X1[1]  = ui().lineEditY1->text().toDouble();
  m_X1[2]  = ui().lineEditZ1->text().toDouble();
  m_X2[0]  = ui().lineEditX2->text().toDouble();
  m_X2[1]  = ui().lineEditY2->text().toDouble();
  m_X2[2]  = ui().lineEditZ2->text().toDouble();
  m_Length = ui().lineEditEL->text().toDouble();
}

double GuiEdgeLengthSourceBox::edgeLength(vec3_t x)
{
  bool outside = false;
  for (int i = 0; i < 3; ++i) {
    if (x[i] < m_X1[i] || x[i] > m_X2[i]) {
      outside = true;
      break;
    }
  }
  if (outside) {
    return -1;
  }
  return m_Length;
}


