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
