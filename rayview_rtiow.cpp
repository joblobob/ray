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
    , cam(QVector3D(0, 0, 1), QVector3D(0, 0, -1), QVector3D(0, 1, 0), (QVector3D(0, 0, 1) - QVector3D(0, 0, -1)).length(), 90, img::aspect_ratio)
{
    ui->setupUi(this);
    m_isColorOnly = false;
    m_shpereY = 0.5f;
    m_testVal2 = 0.0f;

    m_numSamples = 1;
    m_depth = 10;

    scene = new QGraphicsScene(this);

    ui->graphicsView->setScene(scene);

    // materials

    auto material_ground = std::make_shared<lambertian>(QVector3D(0.8, 0.8, 0.0));
    auto material_center = std::make_shared<lambertian>(QVector3D(0.7, 0.3, 0.3));
    auto material_left = std::make_shared<dielectric>(1.5);
    auto material_right = std::make_shared<metal>(QVector3D(0.6, 0.6, 0.6), 0.0);

    // objects

    std::shared_ptr<shape> sph1(new sphere { { 1.0f, 0.0f, -1.0f }, 0.6f, material_right });
    std::shared_ptr<shape> sph3(new sphere { { 0.0f, m_shpereY, -1.0f }, 0.5f, material_center });
    std::shared_ptr<shape> sph4(new sphere { { -1.0f, 0.0f, -1.0f }, 0.6f, material_left });
    std::shared_ptr<shape> sph2(new sphere { { 0.0f, -100.2f, -1.0f }, 100.0f, material_ground });

    auto light = std::make_shared<diffuse_light>(QVector3D(10, 10, 10));

    m_worldObjects.push_back(std::make_shared<xz_rect>(213, 343, 227, 332, 554, light));

    m_worldObjects.push_back(std::move(sph3));
    m_worldObjects.push_back(std::move(sph2));
    m_worldObjects.push_back(std::move(sph1));
    m_worldObjects.push_back(std::move(sph4));

    worldLights = std::make_shared<xz_rect>(213, 343, 227, 332, 554, std::shared_ptr<material>());
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

void RayView_rtiow::renderOneByOne(int width, int height, int samples, const camera& cam, const std::vector<std::shared_ptr<shape>>& worldObjects, int max_depth)
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
                auto u = (row + calc::random_double01()) / (width - 1);
                auto v = (col + calc::random_double01()) / (height - 1);
                Ray r = cam.get_ray(u, v);
                pixel_color += calc::ray_color(r, img::gradientBgVec, worldObjects, worldLights, max_depth, m_isColorOnly);
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
    ui->lblNumSamples->setText("NumSamples: " + QString::number(m_numSamples));
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
