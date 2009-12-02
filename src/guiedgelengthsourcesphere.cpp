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

#include "guiedgelengthsourcesphere.h"

GuiEdgeLengthSourceSphere::GuiEdgeLengthSourceSphere()
{
  m_Centre = vec3_t(0,0,0);
  m_Radius = 1;
  m_Length1 = 1;
  m_Length2 = 1;
}

QString GuiEdgeLengthSourceSphere::write()
{
  QString txt = "sphere: " + name() + "; ";
  QString num;
  num.setNum(m_Centre[0]); txt += num + "; ";
  num.setNum(m_Centre[1]); txt += num + "; ";
  num.setNum(m_Centre[2]); txt += num + "; ";
  num.setNum(m_Radius); txt += num + "; ";
  num.setNum(m_Length1); txt += num + "; ";
  num.setNum(m_Length2); txt += num + "; ";
  return txt;
}

bool GuiEdgeLengthSourceSphere::read(QString txt)
{
  QStringList parts = txt.split(":");
  if (parts.size() != 2) {
    return false;
  }
  if (parts[0].trimmed() == "sphere") {
    QStringList words = parts[1].split(";");
    if (words.size() < 7) {
      return false;
    }
    m_Name = words[0];
    m_Centre = vec3_t(words[1].toDouble(), words[2].toDouble(), words[3].toDouble());
    m_Radius = words[4].toDouble();
    m_Length1 = words[5].toDouble();
    m_Length2 = words[6].toDouble();
  } else {
    return false;
  }
  return true;
}

void GuiEdgeLengthSourceSphere::setDlgFields()
{
  ui().lineEditName->setText(name());
  QString num;
  num.setNum(m_Centre[0]); ui().lineEditX->setText(num);
  num.setNum(m_Centre[1]); ui().lineEditY->setText(num);
  num.setNum(m_Centre[2]); ui().lineEditZ->setText(num);
  num.setNum(m_Radius); ui().lineEditR->setText(num);
  num.setNum(m_Length1); ui().lineEditEL1->setText(num);
  num.setNum(m_Length2); ui().lineEditEL2->setText(num);
}

void GuiEdgeLengthSourceSphere::readDlgFields()
{
  setName(ui().lineEditName->text());
  m_Centre[0] = ui().lineEditX->text().toDouble();
  m_Centre[1] = ui().lineEditY->text().toDouble();
  m_Centre[2] = ui().lineEditZ->text().toDouble();
  m_Radius    = ui().lineEditR->text().toDouble();
  m_Length1   = ui().lineEditEL1->text().toDouble();
  m_Length2   = ui().lineEditEL2->text().toDouble();
}

double GuiEdgeLengthSourceSphere::edgeLength(vec3_t x)
{
  double r = (x - m_Centre).abs();
  if (r > m_Radius) {
    return -1;
  }
  return m_Length1 + r/m_Radius*(m_Length2 - m_Length1);
}

