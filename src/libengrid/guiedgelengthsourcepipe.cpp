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

#include "guiedgelengthsourcepipe.h"

GuiEdgeLengthSourcePipe::GuiEdgeLengthSourcePipe()
{
  m_X1 = vec3_t(0,0,0);
  m_X2 = vec3_t(1,0,0);
  m_R1 = 1;
  m_R2 = 1;
  m_Length1 =1;
  m_Length2 =1;
}

QString GuiEdgeLengthSourcePipe::write()
{
  QString txt = "pipe: " + name() + "; ";
  QString num;
  num.setNum(m_X1[0]); txt += num + "; ";
  num.setNum(m_X1[1]); txt += num + "; ";
  num.setNum(m_X1[2]); txt += num + "; ";
  num.setNum(m_X2[0]); txt += num + "; ";
  num.setNum(m_X2[1]); txt += num + "; ";
  num.setNum(m_X2[2]); txt += num + "; ";
  num.setNum(m_R1); txt += num + "; ";
  num.setNum(m_R2); txt += num + "; ";
  num.setNum(m_Length1); txt += num + "; ";
  num.setNum(m_Length2); txt += num + "; ";
  return txt;
}

bool GuiEdgeLengthSourcePipe::read(QString txt)
{
  QStringList parts = txt.split(":");
  if (parts.size() != 2) {
    return false;
  }
  if (parts[0].trimmed() == "pipe") {
    QStringList words = parts[1].split(";");
    if (words.size() < 11) {
      return false;
    }
    m_Name = words[0].trimmed();
    m_X1 = vec3_t(words[1].toDouble(), words[2].toDouble(), words[3].toDouble());
    m_X2 = vec3_t(words[4].toDouble(), words[5].toDouble(), words[6].toDouble());
    m_R1 = words[7].toDouble();
    m_R2 = words[8].toDouble();
    m_Length1 = words[9].toDouble();
    m_Length2 = words[10].toDouble();
  } else {
    return false;
  }
  return true;
}

void GuiEdgeLengthSourcePipe::setDlgFields()
{
  ui().lineEditName->setText(name());
  QString num;
  num.setNum(m_X1[0]); ui().lineEditX1->setText(num);
  num.setNum(m_X1[1]); ui().lineEditY1->setText(num);
  num.setNum(m_X1[2]); ui().lineEditZ1->setText(num);
  num.setNum(m_X2[0]); ui().lineEditX2->setText(num);
  num.setNum(m_X2[1]); ui().lineEditY2->setText(num);
  num.setNum(m_X2[2]); ui().lineEditZ2->setText(num);
  num.setNum(m_R1); ui().lineEditRi->setText(num);
  num.setNum(m_R2); ui().lineEditRo->setText(num);
  num.setNum(m_Length1); ui().lineEditELi->setText(num);
  num.setNum(m_Length2); ui().lineEditELo->setText(num);
}

void GuiEdgeLengthSourcePipe::readDlgFields()
{
  setName(ui().lineEditName->text());
  m_X1[0]   = ui().lineEditX1->text().toDouble();
  m_X1[1]   = ui().lineEditY1->text().toDouble();
  m_X1[2]   = ui().lineEditZ1->text().toDouble();
  m_X2[0]   = ui().lineEditX2->text().toDouble();
  m_X2[1]   = ui().lineEditY2->text().toDouble();
  m_X2[2]   = ui().lineEditZ2->text().toDouble();
  m_R1      = ui().lineEditRi->text().toDouble();
  m_R2      = ui().lineEditRo->text().toDouble();
  m_Length1 = ui().lineEditELi->text().toDouble();
  m_Length2 = ui().lineEditELo->text().toDouble();
}

double GuiEdgeLengthSourcePipe::edgeLength(vec3_t x)
{
  vec3_t n = m_X2 - m_X1;
  double L = n.abs();
  n.normalise();
  x -= m_X1;
  double l = x*n;
  vec3_t a = l*n;
  double r = (x - a).abs();
  if (r > m_R2 || r < m_R1 || l < 0 || l > L) {
    return -1;
  }
  double w = (r - m_R1)/(m_R2 - m_R1);
  return (1 - w)*m_Length1 + w*m_Length2;
}


