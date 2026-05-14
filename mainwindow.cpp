#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QApplication>
#include <QDebug>

constexpr double mi0 = 4 * M_PI * 1e-7;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //this->setGeometry(300,100,640,480);

    //Initialize canvas object
    m_graphic = new QCustomPlot();
    ui->gridLayout->addWidget(m_graphic, 1, 0, 1, 1);

    ja::HysteresisParams p3C94 = ja::JAHM3::get3C94Params(false);
    m_jamodel = new ja::JAHM3(p3C94, 0.0005, 0.0, {800.0});
    m_jamodel->calculate(ja::SolverType::Euler);

    // Извлечение результатов расчета за один стационарный период
    const auto& B = m_jamodel->getB();
    const auto& H = m_jamodel->getH();

    // Физический анализ полученной петли гистерезиса
    ja::MaterialMetrics metrics = m_jamodel->calculateMetrics();
    double energyLoss = m_jamodel->calculateLosses();
    // Рассчитаем плотность мощности потерь для частоты, например, 50 кГц
    double powerLoss = m_jamodel->calculateTotalPowerLoss(50000.0);

    // Логирование параметров в консоль отладки
    qDebug() << "--- Физические метрики феррита 3C94 ---";
    qDebug() << "Индукция насыщения (Bs):" << metrics.Bs << "T";
    qDebug() << "Коэрцитивная сила (Hc):" << metrics.Hc << "A/m";
    qDebug() << "Остаточная индукция (Br):" << metrics.Br << "T";
    qDebug() << "Энергия потерь за цикл:" << energyLoss << "J/m^3";
    qDebug() << "Мощность потерь на 50 кГц:" << (powerLoss/1000.0) << "kW/m^3";

    QDateTime local(QDateTime::currentDateTime());
    QString filename = QString("bh_data_%1.csv").arg(local.toString("yyyy-MM-dd_hh-mm-ss"));

    auto success = m_jamodel->saveBHToFile(filename);

    if(!success) {
        // In case of an error, an error window is displayed
        QMessageBox::critical(
            nullptr,
            "Failed to save file",
            QString("Failed to save the file:\n%1\n\nPlease check your access permissions or disk space.").arg(filename)
            );
    }

    // Заполнение структуры данных QCustomPlot для отрисовки замкнутой кривой B(H)
    QVector<QCPCurveData> data(H.size());
    m_hyst_curve = new QCPCurve(m_graphic->xAxis, m_graphic->yAxis);
    for(size_t i = 0; i < H.size(); ++i) {
        // Параметры: (индекс_точки, значение_X, значение_Y)
        data[i] = QCPCurveData(static_cast<double>(i), H[i], B[i]);
    }
    m_hyst_curve->setName("B/H Loop");
    m_hyst_curve->data()->set(data, true);
    m_hyst_curve->setPen(QPen(QColor(250, 120, 0), 2)); // Сделаем линию чуть толще (2px)

    // Установка подписей и диапазонов осей
    m_graphic->xAxis->setLabel("Magnetic Field Strength, H (A/m)");
    m_graphic->yAxis->setLabel("Magnetic Flux Density, B (T)");

    double xMin = *std::min_element(H.cbegin(), H.cend());
    double xMax = *std::max_element(H.cbegin(), H.cend());
    double yMin = *std::min_element(B.cbegin(), B.cend());
    double yMax = *std::max_element(B.cbegin(), B.cend());

    m_graphic->xAxis->setRange(xMin, xMax);
    m_graphic->yAxis->setRange(yMin, yMax);
    m_graphic->rescaleAxes();

    // Эстетические отступы по краям холста
    m_graphic->xAxis->scaleRange(1.1, m_graphic->xAxis->range().center());
    m_graphic->yAxis->scaleRange(1.1, m_graphic->yAxis->range().center());

    m_hyst_curve->setVisible(true);

    // Использование QCPItemLine вместо QCPCurve для построения идеальной вертикальной оси H=0
    QCPItemLine *vertical_axis = new QCPItemLine(m_graphic);
    vertical_axis->start->setCoords(0.0, yMin * 1.2); // Начало линии снизу
    vertical_axis->end->setCoords(0.0, yMax * 1.2);   // Конец линии сверху
    vertical_axis->setPen(QPen(Qt::gray, 1, Qt::DashLine));
    vertical_axis->setSelectable(false); // Защита от случайного клика

    // Настройка сетки графика
    m_graphic->xAxis->grid()->setVisible(true);
    m_graphic->yAxis->grid()->setVisible(true);
    m_graphic->xAxis->grid()->setSubGridVisible(true);

    // Включение интерактивности: перетаскивание холста, зум колесиком мыши, выбор объектов
    m_graphic->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    m_graphic->replot();

}

MainWindow::~MainWindow()
{
    delete m_jamodel;
    delete ui;
}

void MainWindow::slotMousePress(QMouseEvent *event)
{
    //Determine the X coordinate on the chart where the mouse click was made
    double coordX = m_graphic->xAxis->pixelToCoord(event->pos().x());
    double coordY = m_graphic->yAxis->pixelToCoord(event->pos().y());

    //qDebug() << coordX;
    //qDebug() << coordY;

    // Берем актуальные границы видимой области графика
    double yMinVisible = m_graphic->yAxis->range().lower;
    double yMaxVisible = m_graphic->yAxis->range().upper;

    QVector<double> x(2), y(2);
    x[0] = coordX;
    y[0] = yMinVisible;
    x[1] = coordX;
    y[1] = yMaxVisible;

    //Set new coordinates
    m_vertical_line->setData(x, y);
    m_vertical_line->setPen(QPen(QColor(250, 120, 0)));


    ui->lineEdit->setText(QString("H: %1 A/m | B: %2 T")
                              .arg(coordX, 0, 'f', 2)
                              .arg(coordY, 0, 'f', 4));
    m_graphic->replot(); //Draw the contents of the canvas
}

void MainWindow::slotMouseMove(QMouseEvent *event)
{
    /** We draw the contents of the canvas,
     *  if when moving the mouse, its button is pressed,
     *  then we call the processing of the mouse coordinates through the click slot
     * */
    if(QApplication::mouseButtons()) slotMousePress(event);
}
