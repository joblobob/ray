#include <QApplication>

#include <fluid.h>
//#include <rayview.h>

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
