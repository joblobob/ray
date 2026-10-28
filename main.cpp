#include <QApplication>
#include <QCoreApplication>

#include <fluid.h>
#include "FluidBenchmark.h"
//#include <rayview.h>

int main(int argc, char* argv[])
{
	for (int i = 1; i < argc; ++i) {
		if (QString::fromLocal8Bit(argv[i]) == QStringLiteral("--benchmark")) {
			QCoreApplication app(argc, argv);
			return runFluidBenchmark(app.arguments());
		}
	}

	QApplication app(argc, argv);

	//RayView* view = new RayView();
	//view->show();

	fluid* view2 = new fluid();
	view2->setupScene();
	view2->show();

	app.exec();
}
