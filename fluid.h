#ifndef FLUID_H
#define FLUID_H



#include "qelapsedtimer.h"

#include <QDebug>
#include <QDialog>
#include <QFile>
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QImage>

import FlipFluid;


namespace Ui {
class fluid;
}

class fluid : public QDialog
{
	Q_OBJECT

public:
	explicit fluid(QWidget* parent = nullptr);

	~fluid() { delete ui; }


	void setupScene();
	void update();


protected:
private slots:
	void on_pauseBtn_clicked()
	{
		m_paused = !m_paused;
		update();
	}
	void on_resetBtn_clicked() { setupScene(); }

	void on_dsBox4_valueChanged(double d)
	{
		//constants::scale    = d;
		//constants::simwidth = (double)constants::maxwidth / constants::scale;
	}


private:
	Ui::fluid* ui;
	QGraphicsScene* m_scene;

	QImage m_imageCanvas;

	QGraphicsPixmapItem* m_sceneItem;
	std::vector<QGraphicsRectItem*> m_gridItems;
	std::vector<QGraphicsEllipseItem*> m_particleItems;
	QGraphicsEllipseItem* m_obstacleItem;

	int m_frameNr;
	bool m_paused;

	FlipFluid m_f;

	void simulate();
	void draw();
	void requestAnimationFrame();



	QElapsedTimer m_timer;
};

#endif