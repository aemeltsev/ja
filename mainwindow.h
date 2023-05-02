#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>
#include <qcustomplot.h>
#include "jahm/jahm.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    // basis: alpha=-0.01; a=12000; Ms=8e5; k=3000; c=0.2;
    double m_alpha = -0.01;
    uint32_t m_a   = 12000;
    uint32_t m_ms  = 8e5;
    uint32_t m_k   = 3000;
    double m_c     = 0.2;

    Ui::MainWindow *ui;
    QCustomPlot    *m_graphic;      // Declare QCustomPlot object
    QCPCurve       *m_vertical_line;     // Declare object for vertical line
    QCPCurve       *m_hyst_curve;
    ja::JAHM       *m_jamodel;

private slots:
    void slotMousePress(QMouseEvent * event);
    void slotMouseMove(QMouseEvent * event);
};

#endif // MAINWINDOW_H
