#include <QApplication>
#include <QDebug>
#include <QFile>
#include <QTextStream>
#include <QVector3D>
#include <iostream>

#include <rayview.h>

int main(int argc, char *argv[]) {

  QApplication app(argc, argv);

  RayView *view = new RayView();
  view->show();
  app.exec();
}
