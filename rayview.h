#ifndef RAYVIEW_H
#define RAYVIEW_H

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
class RayView;
}

class RayView : public QDialog {
    Q_OBJECT

public:
    explicit RayView(QWidget* parent = nullptr);
    ~RayView();

    double hit_object(const QVector3D& centerPoint, const Ray& ray);

    void writeToStream(QTextStream& stream, const QVector3D& pixel, int samples);
    void writeToImg(QImage& img, int x, int y, const QVector3D& pixel, int samples);

    void renderOneByOne(int width, int height, int samples, const camera& cam, const std::vector<std::unique_ptr<shape>>& worldObjects, int max_depth);
    void renderAll(int width, int height, int samples, const camera& cam, const std::vector<std::unique_ptr<shape>>& worldObjects, int max_depth);
    void renderOneRay();

    QVector3D rgbPerSmaples(const QVector3D& pixel, int samples);

private slots:
    void on_dSpin1_valueChanged(double arg1);
    void on_dSpin2_valueChanged(double arg1);
    void on_chkColor_stateChanged(int arg1);
    void on_chkStyle_stateChanged(int arg1);

    void on_horizontalSlider_valueChanged(int value); //shpereY
    void on_horizontalSlider_4_valueChanged(int value);
    void on_horizontalSlider_3_valueChanged(int value);
    void go();

private:
    Ui::RayView* ui;
    QGraphicsScene* scene;

    float m_shpereY;
    float m_testVal2;
    int m_numSamples;
    int m_depth;
    bool m_isColorOnly;
    bool m_isStyleNormal;

    std::vector<std::unique_ptr<shape>> m_worldObjects;

    QImage m_imageCanvas;
    QVector3D m_default_pixel_color;

    QGraphicsPixmapItem* m_sceneItem;

    // camera
    camera cam;
};

#endif // RAYVIEW_H
