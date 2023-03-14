#ifndef FLUID_H
#define FLUID_H

#include <QDebug>
#include <QDialog>
#include <QElapsedTimer>
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QImage>
#include <QPaintEvent>
#include <QTimer>

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
		//	update();
	}
	void on_resetBtn_clicked() { setupScene(); }


private:
	Ui::fluid* ui;
	QGraphicsScene* m_scene;

	QGraphicsPixmapItem* m_sceneItem;
	std::vector<QGraphicsRectItem*> m_gridItems;
	std::vector<QGraphicsEllipseItem*> m_particleItems;
	QGraphicsEllipseItem* m_obstacleItem;

	bool m_paused;

	FlipFluid m_f;

	void simulate();
	void draw();

	int fps = 0;

	QElapsedTimer m_timer;
	QElapsedTimer m_simtimer;
	QElapsedTimer m_setuptimer;
	QTimer m_drawTimer;

	void drawimage();
	QImage m_imageCanvas;
	QGraphicsPixmapItem* m_pixmapSceneItem;
	QGraphicsScene* m_pixmapScene;
	QPainter* m_painter;
	void drawImageToScene();
};

#endif