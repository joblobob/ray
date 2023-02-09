#ifndef FLUID_H
#define FLUID_H

#include <QDebug>
#include <QDialog>
#include <QElapsedTimer>
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>

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


	QElapsedTimer m_timer;
};

#endif