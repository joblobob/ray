#include <QApplication>
#include <QDebug>
#include <QFile>
#include <QTextStream>
#include <QVector3D>
#include <iostream>

#include <rayview.h>
#include <rayview_rtiow.h>

int main(int argc, char* argv[])
{

    QApplication app(argc, argv);

    RayView* view = new RayView();
    view->show();

    RayView_rtiow* view_rtiow = new RayView_rtiow();
    view_rtiow->show();
    app.exec();
}
