#include "rayview.h"

#include <QtConcurrent>
#include <algorithm>

#include "ui_rayview.h"

RayView::RayView(QWidget* parent)
    : QDialog(parent)
    , ui(new Ui::RayView)
    , m_imageCanvas(img::width, img::height, QImage::Format_RGB32)
    , m_default_pixel_color(img::defaultVec)
    , m_sceneItem(nullptr)
{
    ui->setupUi(this);
    m_isColorOnly = false;
    m_shpereY = 0.5f;
    m_testVal2 = 0.0f;
    m_isStyleNormal = true;

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

RayView::~RayView()
{
    delete ui;
}

void RayView::writeToStream(QTextStream& stream, const QVector3D& pixel, int samples)
{
    auto newCol = rgbPerSmaples(pixel, samples);

    // Write the translated [0,255] value of each color component.
    stream << (int)newCol.x() << ' '
           << (int)newCol.y() << ' '
           << (int)newCol.z() << '\n';
}

QVector3D RayView::rgbPerSmaples(const QVector3D& pixel, int samples)
{
    // Divide the color by the number of samples. and gamma-correct for gamma=2.0.
    auto scale = 1.0f / samples;
    auto r = sqrt(scale * pixel.x());
    auto g = sqrt(scale * pixel.y());
    auto b = sqrt(scale * pixel.z());

    return { 255.0f * std::clamp(r, 0.0f, 0.999f), 255.0f * std::clamp(g, 0.0f, 0.999f), 255.0f * std::clamp(b, 0.0f, 0.999f) };
}

void RayView::writeToImg(QImage& img, int x, int y, const QVector3D& pixel, int samples)
{
    if (x < img.width() && y < img.height() && x >= 0 && y >= 0) {
        QVector3D newCol = rgbPerSmaples(pixel, samples);
        // Write the translated [0,255] value of each color component.

        img.setPixelColor(x, y, QColor(newCol.x(), newCol.y(), newCol.z()));
    }
}

void RayView::renderOneByOne(int width, int height, int samples, const camera& cam, const std::vector<std::unique_ptr<shape>>& worldObjects, int max_depth)
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
            RayView::writeToStream(stream, pixel_color, samples);
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
                auto u = (row + calc::random_double()) / (width - 1);
                auto v = (col + calc::random_double()) / (height - 1);
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
    if (m_sceneItem == nullptr) {
        m_sceneItem = scene->addPixmap(QPixmap::fromImage(m_imageCanvas));
        scene->setSceneRect(m_imageCanvas.rect());
    } else
        m_sceneItem->setPixmap(QPixmap::fromImage(m_imageCanvas));

    ui->m_tempsImage->setText("Image: " + QString::number(imageTimer.elapsed()) + " ms");
}

void RayView::renderOneRay()
{
    int randX = calc::random_double(0.0f, img::height + 2);
    int randY = calc::random_double(0.0f, img::width - 1);
    m_default_pixel_color.setX(0.0f);
    m_default_pixel_color.setY(0.0f);
    m_default_pixel_color.setZ(0.0f);
    for (int s = 0; s < m_numSamples; ++s) {
        auto u = (randY + calc::random_double()) / (img::width - 1);
        auto v = (randX + calc::random_double()) / (img::height - 1);
        m_default_pixel_color += calc::ray_color(cam.get_ray(u, v), m_worldObjects, m_depth, m_isColorOnly);
    }
    RayView::writeToImg(m_imageCanvas, img::width - randY - 1, img::height - randX, m_default_pixel_color, m_numSamples);
}

void RayView::go()
{
    QElapsedTimer timer;
    timer.start();
    int rays = 0;
    // render
    if (m_isStyleNormal) {
        renderOneByOne(img::width, img::height, m_numSamples, cam, m_worldObjects, m_depth);
    } else {
        //renderAll(img::width, img::height, m_numSamples, cam, m_worldObjects, m_depth);

        while (timer.elapsed() < 16.666) {
            renderOneRay();
            rays++;
        }
        if (m_sceneItem == nullptr) {
            m_sceneItem = scene->addPixmap(QPixmap::fromImage(m_imageCanvas));
            scene->setSceneRect(m_imageCanvas.rect());
        } else
            m_sceneItem->setPixmap(QPixmap::fromImage(m_imageCanvas));

        QTimer::singleShot(0, this, &RayView::go);
    }
    ui->m_tempsTotal->setText("Total: " + QString::number(timer.elapsed()) + " ms " + QString::number(rays));
}

void RayView::on_dSpin1_valueChanged(double arg1)
{
    m_shpereY = arg1;
    // go();
}

void RayView::on_dSpin2_valueChanged(double arg1)
{
    m_testVal2 = arg1;
    //  go();
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

    //   go();
}

void RayView::on_horizontalSlider_4_valueChanged(int value)
{
    m_numSamples = value;
    //   go();
}

void RayView::on_horizontalSlider_3_valueChanged(int value)
{
    m_depth = value;
    //   go();
}

void RayView::on_chkColor_stateChanged(int arg1)
{
    m_isColorOnly = (bool)arg1;
    //  go();
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
