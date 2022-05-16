#include "rayview.h"

#include <QtConcurrent>
#include <algorithm>

#include "ui_rayview.h"

RayView::RayView(QWidget* parent)
    : QDialog(parent)
    , ui(new Ui::RayView)
    , m_imageCanvas(img::width, img::height, QImage::Format_RGB32)
    , m_default_pixel_color(0.0f, 0.0f, 0.0f)
    , m_bgColor(1.0f, 1.0f, 1.0f)
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

// couleur
QVector3D
RayView::ray_color(const Ray& r,
    const std::vector<std::unique_ptr<shape>>& worldObjects,
    int depth)
{
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0.0f)
        return { 0.0f, 0.0f, 0.0f };

    std::optional<hitPosition> hitRecord = hitFromList(worldObjects, r, 0.001f, calc::infinity);

    if (m_isColorOnly) {
        // juste couleur
        if (hitRecord.has_value()) {
            QVector3D N = (r.at(hitRecord.value().t) - QVector3D(0.0f, 0.0f, -1.0f)).normalized();
            return 0.5f * QVector3D(N.x() + 1.0f, N.y() + 1.0f, N.z() + 1.0f);
        }
    } else {
        // recursive bounce
        if (hitRecord.has_value()) {
            QVector3D target = hitRecord->point + calc::random_in_hemisphere(hitRecord->normal);
            return 0.5f * ray_color({ hitRecord->point, target - hitRecord->point }, worldObjects, depth - 1);
        }
    }

    // pas de hit, background gradient
    QVector3D unit_direction = r.direction.normalized();
    float t = 0.5f * (unit_direction.y() + 1.0f);
    return (1.0f - t) * m_bgColor + t * QVector3D(0.5f, 0.7f, 1.0f);
}

std::optional<hitPosition>
RayView::hitFromList(const std::vector<std::unique_ptr<shape>>& sphereList,
    const Ray& ray,
    double t_min,
    double t_max)
{
    bool hit_anything = false;
    auto closest_so_far = t_max;
    std::optional<hitPosition> retVal;

    for (const auto& object : sphereList) {
        std::optional<hitPosition> didHit = object->hit(ray, t_min, closest_so_far);
        if (didHit.has_value()) {
            hit_anything = true;
            closest_so_far = didHit->t;
            retVal = didHit.value();
        }
    }

    return retVal;
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
                pixel_color += ray_color(r, worldObjects, max_depth);
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
    for (int col = height + 2; col >= 0; --col) {
        for (int row = 0; row < width; ++row) {
            m_default_pixel_color.setX(0.0f);
            m_default_pixel_color.setY(0.0f);
            m_default_pixel_color.setZ(0.0f);
            for (int s = 0; s < samples; ++s) {
                auto u = (row + calc::random_double()) / (width - 1);
                auto v = (col + calc::random_double()) / (height - 1);
                Ray r = cam.get_ray(u, v);
                m_default_pixel_color += ray_color(r, worldObjects, max_depth);
            }

            RayView::writeToImg(m_imageCanvas, width - row - 1, height - col, m_default_pixel_color, samples);
        }
    }

    ui->m_tempsProcess->setText("Process: " + QString::number(processTimer.elapsed()) + " ms");

    QElapsedTimer imageTimer;
    imageTimer.start();
    if (m_sceneItem == nullptr) {
        m_sceneItem = scene->addPixmap(QPixmap::fromImage(m_imageCanvas));
        scene->setSceneRect(m_imageCanvas.rect());
    } else
        m_sceneItem->setPixmap(QPixmap::fromImage(m_imageCanvas));

    ui->m_tempsImage->setText("Image: " + QString::number(imageTimer.elapsed()) + " ms");
}

void RayView::go()
{
    QElapsedTimer timer;
    timer.start();
    // render
    if (m_isStyleNormal) {
        renderOneByOne(img::width, img::height, m_numSamples, cam, m_worldObjects, m_depth);
    } else {
        renderAll(img::width, img::height, m_numSamples, cam, m_worldObjects, m_depth);
    }
    ui->m_tempsTotal->setText("Total: " + QString::number(timer.elapsed()) + " ms");
}

void RayView::on_dSpin1_valueChanged(double arg1)
{
    m_shpereY = arg1;
    go();
}

void RayView::on_dSpin2_valueChanged(double arg1)
{
    m_testVal2 = arg1;
    go();
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

    go();
}

void RayView::on_horizontalSlider_4_valueChanged(int value)
{
    m_numSamples = value;
    go();
}

void RayView::on_horizontalSlider_3_valueChanged(int value)
{
    m_depth = value;
    go();
}

void RayView::on_chkColor_stateChanged(int arg1)
{
    m_isColorOnly = (bool)arg1;
    go();
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
