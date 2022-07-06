#ifndef RayView_rtiow_H
#define RayView_rtiow_H

#include <QDialog>
#include <QFile>
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QImage>
#include <QTextStream>

#include <Ray.h>
#include <calc.h>
#include <camera.h>
#include <hitPosition.h>
#include <shapes.h>

namespace Ui {
class RayView_rtiow;
}

class RayView_rtiow : public QDialog {
    Q_OBJECT

public:
    explicit RayView_rtiow(QWidget* parent = nullptr);
    ~RayView_rtiow();

    double hit_object(const QVector3D& centerPoint, const Ray& ray);

    void writeToStream(QTextStream& stream, const QVector3D& pixel, int samples);

    void renderOneByOne(int width, int height, int samples, const camera& cam, const std::vector<std::shared_ptr<shape>>& worldObjects, int max_depth);

private slots:
    void on_dSpin1_valueChanged(double arg1);
    void on_dSpin2_valueChanged(double arg1);
    void on_chkColor_stateChanged(int arg1);

    void on_horizontalSlider_valueChanged(int value); //shpereY
    void on_horizontalSlider_4_valueChanged(int value);
    void on_horizontalSlider_3_valueChanged(int value);
    void go();

    void on_btnGo_clicked();

private:
    Ui::RayView_rtiow* ui;
    QGraphicsScene* scene;

    float m_shpereY;
    float m_testVal2;
    int m_numSamples;
    int m_depth;
    bool m_isColorOnly;

    std::vector<std::shared_ptr<shape>> m_worldObjects;
    std::shared_ptr<shape> worldLights;

    QImage m_imageCanvas;
    QVector3D m_default_pixel_color;

    QGraphicsPixmapItem* m_sceneItem;

    // camera
    camera cam;
};

#endif // RayView_rtiow_H
