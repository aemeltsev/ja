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

    m_jamodel = new ja::JAHM2();
    m_jamodel->calculate();

    // Getting results
    const auto& B = m_jamodel->getB();
    const auto& H = m_jamodel->getH();
    const auto& M = m_jamodel->getM();
    const auto& t = m_jamodel->getTime();

    //m_jamodel->viewB();
    //m_jamodel->viewH();
    //m_jamodel->viewM();

    QVector<double> BB, HH, MM;
    BB = QVector<double>::fromStdVector(B);
    HH = QVector<double>::fromStdVector(H);
    MM = QVector<double>::fromStdVector(M);

    qDebug() << "Size of HH:" << HH.size();
    qDebug() << "Size of BB:" << BB.size();
    qDebug() << "First 5 elements of HH:" << HH.mid(0, 5);
    qDebug() << "First 5 elements of BB:" << BB.mid(0, 5);

    m_jamodel->saveBHToFile("bh_data.csv", 3);

    m_hyst_curve = new QCPCurve(m_graphic->xAxis, m_graphic->yAxis);
    if (!m_hyst_curve) {
        qDebug() << "m_hyst_curve is null!";
    }

    m_hyst_curve->setName("B/H");
    m_hyst_curve->setData(HH, BB);
    m_hyst_curve->setPen(QPen(QColor(250, 120, 0)));

    //Set labels of the coordinates
    m_graphic->xAxis->setLabel("H");
    m_graphic->yAxis->setLabel("B");

    //Set the maximum and minimum coordinate values
    double xMin = *std::min_element(HH.constBegin(), HH.constEnd());
    double xMax = *std::max_element(HH.constBegin(), HH.constEnd());
    double yMin = *std::min_element(BB.constBegin(), BB.constEnd());
    double yMax = *std::max_element(BB.constBegin(), BB.constEnd());

    m_graphic->xAxis->setRange(xMin, xMax);
    m_graphic->yAxis->setRange(yMin, yMax);

    m_hyst_curve->setVisible(true);  // The curve is visible
    m_hyst_curve->setLayer("main");  // The curve is on the correct layer

    //Initialize a vertical line
    m_vertical_line = new QCPCurve(m_graphic->xAxis, m_graphic->yAxis);

    //Connect the signals from the mouse events on graphic canvas to slots for processing
    connect(m_graphic, &QCustomPlot::mousePress, this, &MainWindow::slotMousePress);
    connect(m_graphic, &QCustomPlot::mouseMove, this, &MainWindow::slotMouseMove);

    //Data vector for vertical line
    QVector<double> x(2) , y(2);
        x[0] = 0;
        y[0] = yMin;
        x[1] = 0;
        y[1] = yMax;

    m_vertical_line->setName("Vertical");       //add name
    m_vertical_line->setData(x, y);             //and set coordinates

    m_graphic->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    //Draw the contents of the canvas
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

    //Preparing the X-coordinates for moving the vertical line
    QVector<double> x(2), y(2);
    x[0] = coordX;
    y[0] = -500000;
    x[1] = coordX;
    y[1] = 500000;

    //Set new coordinates
    m_vertical_line->setData(x, y);
    m_vertical_line->setPen(QPen(QColor(250, 120, 0)));


    ui->lineEdit->setText("H(x): " + QString::number(coordX) +
                          " B(y): " + QString::number(coordY));
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
