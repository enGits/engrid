// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include "guiinflatemesh.h"
#include "guimainwindow.h"

GuiInflateMesh::GuiInflateMesh()
{
}

void GuiInflateMesh::setupSurfaceParameters()
{
  {
    QString buffer;
    QTextStream out(&buffer, QIODevice::WriteOnly);
    out << "\n";
    out << m_EdgeLength << "\n";
    out << "0\n";
    out << "1.5\n";
    out << "0\n";
    out << "1\n";
    out << "1\n";
    GuiMainWindow::pointer()->setXmlSection("engrid/surface/settings", buffer);
  }
  {
    QString buffer;
    QTextStream out(&buffer, QIODevice::WriteOnly);
    out << "\n";
    out << "1 1\n";
    out << "2 " << m_EdgeLength << "\n";
    GuiMainWindow::pointer()->setXmlSection("engrid/surface/table", buffer);
  }
}
