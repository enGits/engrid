#ifndef GUIBRLCADIMPORTDIALOGUE_H
#define GUIBRLCADIMPORTDIALOGUE_H

#include <QDialog>

namespace Ui {
class GuiBrlCadImportDialogue;
}

class GuiBrlCadImportDialogue : public QDialog
{
  Q_OBJECT
  
public:
  explicit GuiBrlCadImportDialogue(QWidget *parent = 0);
  ~GuiBrlCadImportDialogue();
  
private:
  Ui::GuiBrlCadImportDialogue *ui;
};

#endif // GUIBRLCADIMPORTDIALOGUE_H
