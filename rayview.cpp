#include "rayview.h"

#include <QtConcurrent>
#include <algorithm>
#include <execution>

#include <bvh_node.h>
#include <hitPosition.h>

#include <calc.h>

#include "ui_rayview.h"

RayView::RayView(QWidget* parent)
    : QDialog(parent)
    , ui(new Ui::RayView)
    , m_imageCanvas(img::width, img::height, QImage::Format_RGB32)
    , m_default_pixel_color(img::defaultVec)
    , m_lookAt(278, 278, 0)
    , m_lookFrom(278, 278, -800)
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

    //m_worldObjects = random_scene();
    m_worldObjects = box_scene();
    //m_worldObjects = bvh_scene();
    worldLights = std::make_shared<xz_rect>(213, 343, 227, 332, 554, std::shared_ptr<material>());
    //worldLights = std::make_shared<sphere>(QVector3D(190, 90, 190), 90, std::shared_ptr<material>());
}

RayView::~RayView()
{
    delete ui;
}

void RayView::writeToImg(QImage& img, int x, int y, const QVector3D& pixel, int samples)
{
    if (x < img.width() && y < img.height() && x >= 0 && y >= 0) {
        if (m_allPixelsMaps[x][y] < m_numSamples)
            m_allPixelsMaps[x][y]++;
        if (m_allPixelsMaps[x][y] >= m_numSamples)
            m_allPixelsMaps[x][y] = 1;
        QVector3D newCol = calc::rgbPerSamples(pixel, samples);
        // Write the translated [0,255] value of each color component.

        img.setPixelColor(x, y, QColor(newCol.x(), newCol.y(), newCol.z()));
    }
}

void RayView::renderAll(int width, int height, int samples, const camera& cam, const std::vector<std::shared_ptr<shape>>& worldObjects, int max_depth)
{
    QElapsedTimer processTimer;
    processTimer.start();
    bool skip = false;
    int skipped = 0;
    for (int col = height; col >= 0; --col) {
        for (int row = 0; row < width; ++row) {
            skip = false;
            m_default_pixel_color.setX(0.0f);
            m_default_pixel_color.setY(0.0f);
            m_default_pixel_color.setZ(0.0f);
            int x = width - row - 1;
            int y = height - col;
            for (int s = 0; s < samples && !skip; ++s) {
                auto u = (row + calc::random_double01()) / (width - 1);
                auto v = (col + calc::random_double01()) / (height - 1);
                Ray r = cam.get_ray(u, v);
                m_default_pixel_color += calc::ray_color(r, img::gradientBgVec, worldObjects, worldLights, max_depth, m_isColorOnly);
                if (x < m_imageCanvas.width() && y < m_imageCanvas.height() && x >= 0 && y >= 0)
                    if (m_imageCanvas.pixelColor(x, y).red() == m_default_pixel_color.x() && m_imageCanvas.pixelColor(x, y).green() == m_default_pixel_color.y() && m_imageCanvas.pixelColor(x, y).blue() == m_default_pixel_color.z())
                        skip = true;
            }
            if (!skip)
                RayView::writeToImg(m_imageCanvas, x, y, m_default_pixel_color, samples);
            else
                skipped++;
        }
    }

    ui->m_tempsProcess->setText("Process: " + QString::number(processTimer.elapsed()) + " ms");
    qCritical() << skipped;
    QElapsedTimer imageTimer;
    imageTimer.start();

    drawImageToScene();

    ui->m_tempsImage->setText("Image: " + QString::number(imageTimer.elapsed()) + " ms");
}

void RayView::renderOneRay()
{
    int randX = calc::random_double(0.0f, img::width);
    int randY = calc::random_double(0.0f, img::height);

    m_default_pixel_color.setX(0);
    m_default_pixel_color.setY(0);
    m_default_pixel_color.setZ(0);

    int numSamples = 1;
    if (img::width - randX - 1 < m_imageCanvas.width() && img::height - randY - 1 < m_imageCanvas.height() && img::width - randX - 1 >= 0 && img::height - randY - 1 >= 0)
        numSamples = m_allPixelsMaps[img::width - randX - 1][img::height - randY - 1];

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
    QElapsedTimer timer;
    timer.start();
    int rays = 0;
    // render
    if (m_isStyleNormal) {
        renderAll(img::width, img::height, m_numSamples, cam, m_worldObjects, m_depth);
        QTimer::singleShot(timer.elapsed(), this, &RayView::go);
    } else {

        while (timer.elapsed() < 16) {
            renderOneRay();
            rays++;
        }

        drawImageToScene();

        QTimer::singleShot(0, this, &RayView::go);
    }
    ui->m_tempsTotal->setText("Total: " + QString::number(timer.elapsed()) + " ms " + QString::number(rays));
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
        m_lookFrom -= { 0.0f, 0.0f, 5.1f };
        break;
    case Qt::Key_A:
        m_lookFrom += { 5.1f, 0.0f, 0.0f };
        break;
    case Qt::Key_S:
        m_lookFrom += { 0.0f, 0.0f, 5.1f };
        break;
    case Qt::Key_D:
        m_lookFrom -= { 5.1f, 0.0f, 0.0f };
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
    auto light = std::make_shared<diffuse_light>(QVector3D(10, 10, 10));

    world.push_back(std::make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    world.push_back(std::make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    world.push_back(std::make_shared<xz_rect>(213, 343, 227, 332, 554, light));
    world.push_back(std::make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    world.push_back(std::make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    world.push_back(std::make_shared<rectangle>(0, 555, 0, 555, 555, white));

    world.push_back(std::make_shared<box>(QVector3D(130, m_shpereY, 65), QVector3D(295, 165, 230), white));
    world.push_back(std::make_shared<box>(QVector3D(265, m_shpereY, 295), QVector3D(430, 330, 460), white));

    return world;
}

std::vector<std::shared_ptr<shape>> RayView::bvh_scene()
{
    std::vector<std::shared_ptr<shape>> world {};
    static std::vector<std::shared_ptr<shape>> boxes {};
    auto ground = std::make_shared<lambertian>(QVector3D(0.48, 0.83, 0.53));

    const int boxes_per_side = 20;
    for (int i = 0; i < boxes_per_side; i++) {
        for (int j = 0; j < boxes_per_side; j++) {
            auto w = 100.0;
            auto x0 = -1000.0 + i * w;
            auto z0 = -1000.0 + j * w;
            auto y0 = 0.0;
            auto x1 = x0 + w;
            auto y1 = calc::random_double(1, 101);
            auto z1 = z0 + w;

            boxes.push_back(std::make_shared<box>(QVector3D(x0, y0, z0), QVector3D(x1, y1, z1), ground));
        }
    }

    world.push_back(std::make_shared<bvh_node>(boxes, 0, 1));

    auto light = std::make_shared<diffuse_light>(QVector3D(5, 5, 5));
    world.push_back(std::make_shared<xz_rect>(123, 423, 147, 412, 554, light));

    //boule bleu
    auto boundary = std::make_shared<sphere>(QVector3D(360, 150, 145), 70, std::make_shared<dielectric>(1.5));
    world.push_back(boundary);
    world.push_back(std::make_shared<constant_medium>(boundary, 0.2, QVector3D(0.2, 0.4, 0.9)));
    boundary = std::make_shared<sphere>(QVector3D(0, 0, 0), 5000, std::make_shared<dielectric>(1.5));
    world.push_back(std::make_shared<constant_medium>(boundary, .0001, QVector3D(1, 1, 1)));

    // bulle et métal
    world.push_back(std::make_shared<sphere>(QVector3D(260, 150, 45), 50, std::make_shared<dielectric>(1.5)));
    world.push_back(std::make_shared<sphere>(QVector3D(0, 150, 145), 50, std::make_shared<metal>(QVector3D(0.8, 0.8, 0.9), 1.0)));

    //shperes
    std::vector<std::shared_ptr<shape>> boxes2;
    auto white = std::make_shared<lambertian>(QVector3D(.73, .73, .73));
    int ns = 100;
    for (int j = 0; j < ns; j++) {
        boxes2.push_back(std::make_shared<sphere>(QVector3D(calc::randomVec(1, 165)), 10, white));
    }

    world.push_back(std::make_shared<bvh_node>(boxes2, 0, 1));

    return world;
}
