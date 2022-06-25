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
    , cam(m_lookFrom, m_lookAt, m_vup)
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

    // materials

    auto material_ground = std::make_shared<lambertian>(QVector3D(0.8, 0.8, 0.0));
    auto material_center = std::make_shared<lambertian>(QVector3D(0.7, 0.3, 0.3));
    // auto material_left = std::make_shared<metal>(QVector3D(0.8, 0.8, 0.8), 0.3);
    //  auto material_right = std::make_shared<metal>(QVector3D(0.8, 0.6, 0.2), 1.0);
    auto material_left = std::make_shared<dielectric>(1.5);
    auto material_right = std::make_shared<metal>(QVector3D(0.8, 0.6, 0.2), 0.0);

    // objects

    std::unique_ptr<shape> sph1(new sphere { { 1.0f, 0.0f, -1.0f }, 0.6f, material_right });
    std::unique_ptr<shape> sph3(new sphere { { 0.0f, m_shpereY, -1.0f }, 0.5f, material_center });
    std::unique_ptr<shape> sph4(new sphere { { -1.0f, 0.0f, -1.0f }, 0.6f, material_left });
    std::unique_ptr<shape> sph2(new sphere { { 0.0f, -100.2f, -1.0f }, 100.0f, material_ground });

    // std::unique_ptr<shape> rect1(new rectangle { { 0.0f, m_testVal, -20.0f }, 1.0f, 1.0f, 1.0f });
    m_worldObjects.push_back(std::move(sph3));
    m_worldObjects.push_back(std::move(sph2));
    m_worldObjects.push_back(std::move(sph1));
    m_worldObjects.push_back(std::move(sph4));
    //worldObjects.push_back(std::move(rect1));
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

void RayView::renderAll(int width, int height, int samples, const camera& cam, const std::vector<std::unique_ptr<shape>>& worldObjects, int max_depth)
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
                m_default_pixel_color += calc::ray_color(r, worldObjects, max_depth, m_isColorOnly);
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

        m_default_pixel_color += calc::ray_color(cam.get_ray(u, v), m_worldObjects, m_depth, m_isColorOnly);
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
        m_lookFrom -= { 0.0f, 0.0f, 0.1f };
        break;
    case Qt::Key_A:
        m_lookFrom += { 0.1f, 0.0f, 0.0f };
        break;
    case Qt::Key_S:
        m_lookFrom += { 0.0f, 0.0f, 0.1f };
        break;
    case Qt::Key_D:
        m_lookFrom -= { 0.1f, 0.0f, 0.0f };
        break;
    }

    cam.resetCam(m_lookFrom, m_lookAt, m_vup);
}
