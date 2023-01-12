#include <QApplication>
#include <QDebug>
#include <QFile>
#include <QTextStream>
#include <QVector3D>
#include <iostream>

#include <rayview.h>
#include <fluid.h>

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

    //RayView* view = new RayView();
    //view->show();

    fluid* view2 = new fluid();
	view2->setupScene();
	view2->show();

    app.exec();
}
