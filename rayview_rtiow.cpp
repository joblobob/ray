#include "rayview_rtiow.h"

#include <QtConcurrent>
#include <algorithm>

#include "ui_rayview_rtiow.h"

RayView_rtiow::RayView_rtiow(QWidget* parent)
    : QDialog(parent)
    , ui(new Ui::RayView_rtiow)
    , m_imageCanvas(img::width, img::height, QImage::Format_RGB32)
    , m_default_pixel_color(img::defaultVec)
    , m_sceneItem(nullptr)
{
    ui->setupUi(this);
    m_isColorOnly = false;
    m_shpereY = 0.5f;
    m_testVal2 = 0.0f;

    m_numSamples = 1;
    m_depth = 10;

    scene = new QGraphicsScene(this);

    ui->graphicsView->setScene(scene);

    // objects

    std::unique_ptr<shape> sph1(new sphere { { 0.0f, m_shpereY, -1.0f }, 0.6f });
    std::unique_ptr<shape> sph3(new sphere { { 0.0f, m_shpereY, -1.5f }, 0.1f });
    std::unique_ptr<shape> sph2(new sphere { { 0.0f, -100.2f, -1.0f }, 100.0f });
    // std::unique_ptr<shape> rect1(new rectangle { { 0.0f, m_testVal, -20.0f }, 1.0f, 1.0f, 1.0f });
    m_worldObjects.push_back(std::move(sph1));
    m_worldObjects.push_back(std::move(sph2));
    m_worldObjects.push_back(std::move(sph3));
    //worldObjects.push_back(std::move(rect1));
}

RayView_rtiow::~RayView_rtiow()
{
    delete ui;
}

void RayView_rtiow::writeToStream(QTextStream& stream, const QVector3D& pixel, int samples)
{
    auto newCol = calc::rgbPerSamples(pixel, samples);

    // Write the translated [0,255] value of each color component.
    stream << (int)newCol.x() << ' '
           << (int)newCol.y() << ' '
           << (int)newCol.z() << '\n';
}

void RayView_rtiow::renderOneByOne(int width, int height, int samples, const camera& cam, const std::vector<std::unique_ptr<shape>>& worldObjects, int max_depth)
{
    QElapsedTimer processTimer;
    processTimer.start();

    QByteArray data;
    QTextStream stream(&data);
    stream << "P3\n"
           << width << ' ' << height << "\n255\n";
    QVector3D pixel_color(0, 0, 0);
    for (int col = height + 2; col >= 0; --col) {
        for (int row = 0; row < width; ++row) {
            pixel_color.setX(0);
            pixel_color.setY(0);
            pixel_color.setZ(0);
            for (int s = 0; s < samples; ++s) {
                auto u = (row + calc::random_double()) / (width - 1);
                auto v = (col + calc::random_double()) / (height - 1);
                Ray r = cam.get_ray(u, v);
                pixel_color += calc::ray_color(r, worldObjects, max_depth, m_isColorOnly);
            }
            RayView_rtiow::writeToStream(stream, pixel_color, samples);
        }
    }

    ui->m_tempsProcess->setText("Process: " + QString::number(processTimer.elapsed()) + " ms");

    QElapsedTimer imageTimer;
    imageTimer.start();
    QPixmap image;
    image.loadFromData(data);
    scene->addPixmap(image);
    scene->setSceneRect(image.rect());
    ui->m_tempsImage->setText("Image: " + QString::number(imageTimer.elapsed()) + " ms");
}

void RayView_rtiow::go()
{
    QElapsedTimer timer;
    timer.start();

    // render
    renderOneByOne(img::width, img::height, m_numSamples, cam, m_worldObjects, m_depth);

    ui->m_tempsTotal->setText("Total: " + QString::number(timer.elapsed()) + " ms " + QString::number(img::width * img::height));
}

void RayView_rtiow::on_dSpin1_valueChanged(double arg1)
{
    m_shpereY = arg1;
}

void RayView_rtiow::on_dSpin2_valueChanged(double arg1)
{
    m_testVal2 = arg1;
}

//sphereY
void RayView_rtiow::on_horizontalSlider_valueChanged(int value)
{
    m_shpereY = value * 0.01f;
    if (!m_worldObjects.empty())
        m_worldObjects.at(0)->center.setY(m_shpereY);

    ui->dSpin1->blockSignals(true);
    ui->dSpin1->setValue(m_shpereY);
    ui->dSpin1->blockSignals(false);
}

void RayView_rtiow::on_horizontalSlider_4_valueChanged(int value)
{
    m_numSamples = value;
}

void RayView_rtiow::on_horizontalSlider_3_valueChanged(int value)
{
    m_depth = value;
}

void RayView_rtiow::on_chkColor_stateChanged(int arg1)
{
    m_isColorOnly = (bool)arg1;
}

void RayView_rtiow::on_btnGo_clicked()
{
    go();
}
