#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QApplication>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //this->setGeometry(300,100,640,480);

    //Initialize canvas object
    m_graphic = new QCustomPlot();
    ui->gridLayout->addWidget(m_graphic, 1, 0, 1, 1);

    //Initialize a vertical line
    m_vertical_line = new QCPCurve(m_graphic->xAxis, m_graphic->yAxis);

    //Connect the signals from the mouse events on graphic canvas to slots for processing
    connect(m_graphic, &QCustomPlot::mousePress, this, &MainWindow::slotMousePress);
    connect(m_graphic, &QCustomPlot::mouseMove, this, &MainWindow::slotMouseMove);

    //Data vector for vertical line
    QVector<double> x(2) , y(2);
        x[0] = 0;
        y[0] = -500000;
        x[1] = 0;
        y[1] = 500000;

    //m_graphic->addPlottable(m_vertical_line);   //Add line to canvas
    m_vertical_line->setName("Vertical");       //add name
    m_vertical_line->setData(x, y);             //and set coordinates

    m_jamodel = new ja::JAHM(m_ms, m_a, m_k, m_alpha, m_c);
    m_jamodel->jaTotalMagnetizeCalc(0, 40, 0.004);
    QVector<double> h , m;
    h = QVector<double>::fromStdVector(m_jamodel->jaGetH());
    m = QVector<double>::fromStdVector(m_jamodel->jaGetM());

    m_hyst_curve = new QCPCurve(m_graphic->xAxis, m_graphic->yAxis);
    m_hyst_curve->setName("M/H");
    m_hyst_curve->setData(h, m);

    //Set labels of the coordinates
    m_graphic->xAxis->setLabel("x");
    m_graphic->yAxis->setLabel("y");

    //Set the maximum and minimum coordinate values
    m_graphic->xAxis->setRange(-500000, 500000);
    m_graphic->yAxis->setRange(-500000, 500000);

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

    qDebug() << coordX;
    qDebug() << coordY;

    //Preparing the X-coordinates for moving the vertical line
    QVector<double> x(2), y(2);
    x[0] = coordX;
    y[0] = -500000;
    x[1] = coordX;
    y[1] = 500000;

    //Set new coordinates
    m_vertical_line->setData(x, y);
    m_vertical_line->setPen(QPen(QColor(250, 120, 0)));


    ui->lineEdit->setText("x: " + QString::number(coordX) +
                          " y: " + QString::number(coordY));
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
