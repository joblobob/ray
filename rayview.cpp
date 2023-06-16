#include "rayview.h"

#include <QtConcurrent>
#include <algorithm>
#include <execution>

#include <hitPosition.h>

#include <calc.h>

#include "ui_rayview.h"

RayView::RayView(QWidget* parent)
    : QDialog(parent)
    , ui(new Ui::RayView)
    , m_imageCanvas(img::width, img::height, QImage::Format_RGB32)
    , m_default_pixel_color(img::defaultVec)
    , m_lookAt(0.0f, 0.0f, -1.0f)
    , m_lookFrom(0.0f, 0.0f, 1.0f)
    , m_vup(0.0f, 1.0f, 0.0f)
    , m_sceneItem(nullptr)
    , m_nextFrameX(0)
    , m_nextFrameY(0)
    , m_raysSinceEpoch(0)
    , m_allisDone(false)
    , m_allPixelsMaps { img::width + 1, QVector<int>(img::height + 1, 1) }
    , cam(m_lookFrom, m_lookAt, m_vup, (m_lookFrom - m_lookAt).length())
{
    ui->setupUi(this);
    m_isColorOnly = false;
    m_shpereY = 0.5f;
    m_testVal2 = 0.0f;
    m_isStyleNormal = true;
    m_imageCanvas.fill(Qt::black);

    m_allPixels.reserve(img::width * img::height);

    m_numSamples = 1;
    m_depth = 10;

    scene = new QGraphicsScene(this);

    ui->graphicsView->setScene(scene);

    m_worldObjects = normal_scene();
    //m_worldObjects = box_scene();
	//m_worldObjects = random_scene();
    worldLights = std::make_shared<sphere>(QVector3D(275, 2000, 1000), 500, std::shared_ptr<material>());
    worldLights = std::make_shared<sphere>(QVector3D(190, 90, 190), 90, std::shared_ptr<material>());
}

RayView::~RayView()
{
    delete ui;
}

void RayView::writeToImg(QImage& img, int x, int y, const QVector3D& pixel, int samples)
{
    if (x < img.width() && y < img.height() && x >= 0 && y >= 0) {
        QVector3D newCol = calc::rgbPerSamples(pixel, samples);
        // Write the translated [0,255] value of each color component.

        img.setPixelColor(x, y, QColor(newCol.x(), newCol.y(), newCol.z()));
    }
}

void RayView::renderOneRay()
{
    int randX = calc::random_doubleW();
    int randY = calc::random_doubleH();

    m_default_pixel_color.setX(0);
    m_default_pixel_color.setY(0);
    m_default_pixel_color.setZ(0);

    int numSamples = m_numSamples;

    //random mais dans une liste detat pour faire une passe a 1, ensuite une passe a 2 et faire chaque pixels mais randomly wow
    for (int s = 0; s < numSamples; ++s) {
        float u = (randX + calc::random_double01()) / img::width;
        float v = (randY + calc::random_double01()) / img::height;

        m_default_pixel_color += calc::ray_color(cam.get_ray(u, v), img::gradientBgVec, m_worldObjects, worldLights, m_depth, m_isColorOnly);
    }

    RayView::writeToImg(m_imageCanvas, img::width - randX - 1, img::height - randY - 1, m_default_pixel_color, numSamples);
}

void RayView::drawImageToScene()
{
    if (m_sceneItem == nullptr) {
        m_sceneItem = scene->addPixmap(QPixmap::fromImage(m_imageCanvas));
        scene->setSceneRect(m_imageCanvas.rect());
    } else
        m_sceneItem->setPixmap(QPixmap::fromImage(m_imageCanvas));
}

void RayView::go()
{
    timer.start(); 
    int rays = 0;

    // render
    while (timer.elapsed() < 16) {
        renderOneRay();
        rays++;
    }

    //draw
    drawImageToScene();

	ui->m_tempsTotal->setText("Total: " + QString::number(timer.elapsed()) + " ms " + QString::number(rays));

    //callback in 0ms
    QTimer::singleShot(0, this, &RayView::go);
    
}

void RayView::on_dSpin1_valueChanged(double arg1)
{
    m_shpereY = arg1;
}

void RayView::on_dSpin2_valueChanged(double arg1)
{
    m_testVal2 = arg1;
}

//sphereY
void RayView::on_horizontalSlider_valueChanged(int value)
{
    m_shpereY = value * 0.01f;
    if (!m_worldObjects.empty())
        m_worldObjects.at(0)->center.setY(m_shpereY);

    ui->dSpin1->blockSignals(true);
    ui->dSpin1->setValue(m_shpereY);
    ui->dSpin1->blockSignals(false);
}

void RayView::on_horizontalSlider_4_valueChanged(int value)
{
    m_numSamples = value;
    ui->lblNumSamples->setText("NumSamples: " + QString::number(m_numSamples));
}

void RayView::on_horizontalSlider_3_valueChanged(int value)
{
    m_depth = value;
}

void RayView::on_chkColor_stateChanged(int arg1)
{
    m_isColorOnly = (bool)arg1;
}
void RayView::on_chkStyle_stateChanged(int arg1)
{
    m_isStyleNormal = (bool)arg1;
    if (m_sceneItem != nullptr) {
        delete m_sceneItem;
        m_sceneItem = nullptr;
    }
    go();
}

void RayView::keyPressEvent(QKeyEvent* keyEvent)
{
    switch (keyEvent->key()) {
    case Qt::Key_W:
        m_lookFrom -= { 0.0f, 0.0f, 1.0f };
        break;
    case Qt::Key_A:
        m_lookFrom += { 1.0f, 0.0f, 0.0f };
        break;
    case Qt::Key_S:
        m_lookFrom += { 0.0f, 0.0f, 1.0f };
        break;
    case Qt::Key_D:
        m_lookFrom -= { 1.0f, 0.0f, 0.0f };
        break;
    }

    cam.resetCam(m_lookFrom, m_lookAt, m_vup, (m_lookFrom - m_lookAt).length());
}

std::vector<std::shared_ptr<shape>> RayView::random_scene()
{
    std::vector<std::shared_ptr<shape>> world;

    auto ground_material = std::make_shared<lambertian>(QVector3D(0.5, 0.5, 0.5));
    world.push_back(std::make_shared<sphere>(QVector3D(0, -1000, 0), 1000, ground_material));

    for (int a = -5; a < 5; a++) {
        for (int b = -5; b < 5; b++) {
            auto choose_mat = calc::random_double01();
            QVector3D center(a + 0.9f * calc::random_double01(), 0.2f, b + 0.9f * calc::random_double01());

            if ((center - QVector3D(4.0f, 0.2f, 0.0f)).length() > 0.9f) {
                std::shared_ptr<material> sphere_material;

                if (choose_mat < 0.8f) {
                    // diffuse
                    auto albedo = calc::randomVec() * calc::randomVec();
                    sphere_material = std::make_shared<lambertian>(albedo);
                    world.push_back(std::make_shared<sphere>(center, 0.2f, sphere_material));
                } else if (choose_mat < 0.95f) {
                    // metal
                    auto albedo = calc::randomVec(0.5, 1);
                    auto fuzz = calc::random_double(0, 0.5);
                    sphere_material = std::make_shared<metal>(albedo, fuzz);
                    world.push_back(std::make_shared<sphere>(center, 0.2f, sphere_material));
                } else {
                    // glass
                    sphere_material = std::make_shared<dielectric>(1.5);
                    world.push_back(std::make_shared<sphere>(center, 0.2f, sphere_material));
                }
            }
        }
    }

    auto material1 = std::make_shared<dielectric>(1.5);
    world.push_back(std::make_shared<sphere>(QVector3D(0, 1, 0), 1.0, material1));

    auto material2 = std::make_shared<lambertian>(QVector3D(0.4, 0.2, 0.1));
    world.push_back(std::make_shared<sphere>(QVector3D(-4, 1, 0), 1.0, material2));

    auto material3 = std::make_shared<metal>(QVector3D(0.7, 0.6, 0.5), 0.0);
    world.push_back(std::make_shared<sphere>(QVector3D(4, 1, 0), 1.0, material3));

    return world;
}

std::vector<std::shared_ptr<shape>> RayView::normal_scene()
{
    std::vector<std::shared_ptr<shape>> world;

    // materials

    auto material_ground = std::make_shared<lambertian>(QVector3D(0.8, 0.8, 0.0));
    auto material_center = std::make_shared<lambertian>(QVector3D(0.7, 0.3, 0.3));
    auto material_left = std::make_shared<dielectric>(1.5);
    auto material_right = std::make_shared<metal>(QVector3D(0.6, 0.6, 0.6), 0.0);
    auto difflight = std::make_shared<diffuse_light>(QVector3D(2, 2, 2));

    // objects

    std::shared_ptr<shape> sph1(std::make_shared<sphere>(QVector3D { 1.0f, 0.0f, -1.0f }, 0.6f, material_right));
    std::shared_ptr<shape> sph3(std::make_shared<sphere>(QVector3D { 0.0f, m_shpereY, -1.0f }, 0.5f, material_center));
    std::shared_ptr<shape> sph4(std::make_shared<sphere>(QVector3D { -1.0f, 0.0f, -1.0f }, 0.6f, material_left));
    std::shared_ptr<shape> sph2(std::make_shared<sphere>(QVector3D { 0.0f, -100.2f, -1.0f }, 100.0f, material_ground));

    world.push_back(std::make_shared<sphere>(QVector3D { 0.0f, m_shpereY, -1.0f }, 0.5f, material_center));

    world.push_back(std::move(sph2));
    world.push_back(std::move(sph1));
    world.push_back(std::move(sph4));

    world.push_back(std::make_shared<rectangle>(3, 5, 1, 3, -2, difflight));
    return world;
}

std::vector<std::shared_ptr<shape>> RayView::box_scene()
{
    std::vector<std::shared_ptr<shape>> world;

    auto red = std::make_shared<lambertian>(QVector3D(.65, .05, .05));
    auto white = std::make_shared<lambertian>(QVector3D(.73, .73, .73));
    auto green = std::make_shared<lambertian>(QVector3D(.12, .45, .15));
    auto gray = std::make_shared<lambertian>(QVector3D(.77, .77, .77));
    auto light = std::make_shared<diffuse_light>(QVector3D(20, 20, 20));
    auto metallike = std::make_shared<metal>(QVector3D(0.6, 0.6, 0.6), 0.0);

    world.push_back(std::make_shared<yz_rect>(0, 600, -600, 3000, 600, red));
    world.push_back(std::make_shared<yz_rect>(0, 600, -600, 3000, 0, green));

    world.push_back(std::make_shared<sphere>(QVector3D(275, 3000, 1000), 500, light));

    world.push_back(std::make_shared<xz_rect>(0, 600, -600, 3000, 0, metallike));
    world.push_back(std::make_shared<xz_rect>(0, 200, -600, 3000, 600, white));
    world.push_back(std::make_shared<xz_rect>(400, 600, -600, 3000, 600, white));

    world.push_back(std::make_shared<xz_rect>(200, 400, -600, -200, 580, gray));
    world.push_back(std::make_shared<xz_rect>(200, 400, 200, 400, 580, gray));
    world.push_back(std::make_shared<xz_rect>(200, 400, 1000, 1200, 580, gray));
    world.push_back(std::make_shared<xz_rect>(200, 400, 1800, 2000, 580, gray));
    world.push_back(std::make_shared<xz_rect>(200, 400, 2600, 2800, 580, gray));

    world.push_back(std::make_shared<rectangle>(-600, 600, 0, 600, 3000, white));
    world.push_back(std::make_shared<rectangle>(-600, 600, 0, 600, -1000, light));

    // cpp north
    QVector3D initialPosition = QVector3D(0, 250, 700);

    float scaleFactor = 0.22f;

    //C
    world.push_back(std::make_shared<box>(initialPosition, initialPosition + (QVector3D(100, 200, 100) *= scaleFactor), red));
    world.push_back(
        std::make_shared<box>(initialPosition + (QVector3D(100, 0, 0) *= scaleFactor), initialPosition + (QVector3D(200, 50, 100) *= scaleFactor), red));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 150, 0) *= scaleFactor), initialPosition + (QVector3D(200, 200, 100) *= scaleFactor), red));

    initialPosition += QVector3D(300, 0, 0) *= scaleFactor;

    //P
    world.push_back(
        std::make_shared<box>(initialPosition + (QVector3D(0, 0, 0) *= scaleFactor), initialPosition + (QVector3D(100, 200, 100) *= scaleFactor), red));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 150, 0) *= scaleFactor), initialPosition + (QVector3D(200, 200, 100) *= scaleFactor), red));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 50, 0) *= scaleFactor), initialPosition + (QVector3D(200, 100, 100) *= scaleFactor), red));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(150, 100, 0) *= scaleFactor), initialPosition + (QVector3D(200, 150, 100) *= scaleFactor), red));

    initialPosition += QVector3D(300, 0, 0) *= scaleFactor;

    //P
    world.push_back(
        std::make_shared<box>(initialPosition + (QVector3D(0, 0, 0) *= scaleFactor), initialPosition + (QVector3D(100, 200, 100) *= scaleFactor), red));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 150, 0) *= scaleFactor), initialPosition + (QVector3D(200, 200, 100) *= scaleFactor), red));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 50, 0) *= scaleFactor), initialPosition + (QVector3D(200, 100, 100) *= scaleFactor), red));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(150, 100, 0) *= scaleFactor), initialPosition + (QVector3D(200, 150, 100) *= scaleFactor), red));

    initialPosition += QVector3D(350, 0, 0) *= scaleFactor;

    //N
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(0, 0, 0) *= scaleFactor), initialPosition + (QVector3D(100, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 125, 0) *= scaleFactor), initialPosition + (QVector3D(250, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(150, 0, 0) *= scaleFactor), initialPosition + (QVector3D(250, 150, 100) *= scaleFactor), green));

    initialPosition += QVector3D(350, 0, 0) *= scaleFactor;

    //O
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(0, 0, 0) *= scaleFactor), initialPosition + (QVector3D(100, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 150, 0) *= scaleFactor), initialPosition + (QVector3D(300, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(200, 0, 0) *= scaleFactor), initialPosition + (QVector3D(300, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 0, 0) *= scaleFactor), initialPosition + (QVector3D(300, 50, 100) *= scaleFactor), green));

    initialPosition += QVector3D(400, 0, 0) *= scaleFactor;

    //R
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(0, 0, 0) *= scaleFactor), initialPosition + (QVector3D(100, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 150, 0) *= scaleFactor), initialPosition + (QVector3D(200, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(150, 100, 0) *= scaleFactor), initialPosition + (QVector3D(200, 150, 100) *= scaleFactor), green));

    initialPosition += QVector3D(300, 0, 0) *= scaleFactor;

    //T
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(75, 0, 0) *= scaleFactor), initialPosition + (QVector3D(125, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(0, 150, 0) *= scaleFactor), initialPosition + (QVector3D(200, 200, 100) *= scaleFactor), green));

    initialPosition += QVector3D(300, 0, 0) *= scaleFactor;

    //H
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(0, 0, 0) *= scaleFactor), initialPosition + (QVector3D(100, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(200, 0, 0) *= scaleFactor), initialPosition + (QVector3D(300, 200, 100) *= scaleFactor), green));
    world.push_back(std::make_shared<box>(
        initialPosition + (QVector3D(100, 75, 0) *= scaleFactor), initialPosition + (QVector3D(200, 125, 100) *= scaleFactor), green));

    return world;
}
