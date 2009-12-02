#include <QApplication>
#include <QPushButton>

#include <QHBoxLayout>
#include <QFormLayout>
#include <QtDebug>
#include <QMainWindow>

#include "dialoglineedit.h"

int main(int argc, char *argv[])
{
  QApplication app(argc, argv);

//     QMainWindow window;
  QWidget window;

  QFormLayout* layout = new QFormLayout();
  /*    layout->setSpacing(-1);
      this->setContentsMargins(0,0,0,0);*/
  window.setLayout(layout);

/////////////////////
//my widget
/////////////////////

  DialogLineEdit* hello = new DialogLineEdit();
  hello->setMode(true);
  hello->setCaption("YOYOYOOYOYOYY");
  layout->addWidget(hello);

  DialogLineEdit* hello2 = new DialogLineEdit();
  hello2->setMode(false);
  hello2->setCaption("YOYOYOOYOYOYY");
  layout->addWidget(hello2);

/////////////////////
// a normal line edit
/////////////////////
  /*
      QLineEdit* hello = new QLineEdit();

      layout->addWidget(hello);

      QLineEdit* hello2 = new QLineEdit();
      layout->addWidget(hello2);
  */
/////////////////////

//     DialogLineEdit helloo(&window);
//     window.addWidget(hello);

//      QPushButton quit("Quit", &window);
//      quit.setFont(QFont("Times", 18, QFont::Bold));
//      helloo.setGeometry(10, 40, 180, 40);
//      QObject::connect(&quit, SIGNAL(clicked()), &app, SLOT(quit()));

  window.show();
  return app.exec();
}
