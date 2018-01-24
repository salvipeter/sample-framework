#include <QtWidgets/QApplication>

#include "MyWindow.h"

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  MyWindow window(&app);
  window.show();
  return app.exec();
}
