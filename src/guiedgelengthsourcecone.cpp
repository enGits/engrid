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

#include "guiedgelengthsourcecone.h"

GuiEdgeLengthSourceCone::GuiEdgeLengthSourceCone()
{
  m_X1 = vec3_t(0,0,0);
  m_X2 = vec3_t(1,0,0);
  m_R1 = 1;
  m_R2 = 1;
  m_Length1 =1;
  m_Length2 =1;
}

QString GuiEdgeLengthSourceCone::write()
{
  QString txt = "cone: " + name() + "; ";
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

bool GuiEdgeLengthSourceCone::read(QString txt)
{
  QStringList parts = txt.split(":");
  if (parts.size() != 2) {
    return false;
  }
  if (parts[0].trimmed() == "cone") {
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

void GuiEdgeLengthSourceCone::setDlgFields()
{
  ui().lineEditName->setText(name());
  QString num;
  num.setNum(m_X1[0]); ui().lineEditX1->setText(num);
  num.setNum(m_X1[1]); ui().lineEditY1->setText(num);
  num.setNum(m_X1[2]); ui().lineEditZ1->setText(num);
  num.setNum(m_X2[0]); ui().lineEditX2->setText(num);
  num.setNum(m_X2[1]); ui().lineEditY2->setText(num);
  num.setNum(m_X2[2]); ui().lineEditZ2->setText(num);
  num.setNum(m_R1); ui().lineEditR1->setText(num);
  num.setNum(m_R2); ui().lineEditR2->setText(num);
  num.setNum(m_Length1); ui().lineEditEL1->setText(num);
  num.setNum(m_Length2); ui().lineEditEL2->setText(num);
}

void GuiEdgeLengthSourceCone::readDlgFields()
{
  setName(ui().lineEditName->text());
  m_X1[0]   = ui().lineEditX1->text().toDouble();
  m_X1[1]   = ui().lineEditY1->text().toDouble();
  m_X1[2]   = ui().lineEditZ1->text().toDouble();
  m_X2[0]   = ui().lineEditX2->text().toDouble();
  m_X2[1]   = ui().lineEditY2->text().toDouble();
  m_X2[2]   = ui().lineEditZ2->text().toDouble();
  m_R1      = ui().lineEditR1->text().toDouble();
  m_R2      = ui().lineEditR2->text().toDouble();
  m_Length1 = ui().lineEditEL1->text().toDouble();
  m_Length2 = ui().lineEditEL2->text().toDouble();
}

double GuiEdgeLengthSourceCone::edgeLength(vec3_t x)
{
  vec3_t n = m_X2 - m_X1;
  double L = n.abs();
  n.normalise();
  x -= m_X1;
  double l = x*n;
  double R = m_R1 + l/L*(m_R2 - m_R1);
  vec3_t a = l*n;
  double r = (x - a).abs();
  if (r > R || l < 0 || l > L) {
    return -1;
  }
  return m_Length1 + r/R*(m_Length2 - m_Length1);
}


