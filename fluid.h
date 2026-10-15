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
#include <QResizeEvent>
#include <QTimer>
#include <memory>


namespace Ui {
class fluid;
}

class fluid : public QDialog
{
	Q_OBJECT

public:
	explicit fluid(QWidget* parent = nullptr);

	~fluid();


	void setupScene();
	void update();


protected:
	void resizeEvent(QResizeEvent* event) override;

private slots:
	void on_pauseBtn_clicked()
	{
		m_paused = !m_paused;
		//	update();
	}
	void on_resetBtn_clicked() { setupScene(); }


private:
	Ui::fluid* ui;
	QGraphicsEllipseItem* m_obstacleItem;

	bool m_paused;

	struct Impl;
	std::unique_ptr<Impl> m_impl;

	void simulate();
	void draw();

	int fps = 0;

	QElapsedTimer m_timer;
	QElapsedTimer m_simtimer;
	QElapsedTimer m_setuptimer;
	QTimer m_drawTimer;

	void drawimage();
	QPixmap m_pointsPixmap;
	QGraphicsPixmapItem* m_pixmapSceneItem;
	QGraphicsScene* m_pixmapScene;
	QPainter* m_painter;
	void drawImageToScene();
};

#endif
