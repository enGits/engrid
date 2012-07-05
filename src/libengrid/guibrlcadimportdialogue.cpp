#include "guibrlcadimportdialogue.h"
#include "ui_guibrlcadimportdialogue.h"

GuiBrlCadImportDialogue::GuiBrlCadImportDialogue(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::GuiBrlCadImportDialogue)
{
  ui->setupUi(this);
}

GuiBrlCadImportDialogue::~GuiBrlCadImportDialogue()
{
  delete ui;
}
