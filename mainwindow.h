#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>
#include <qcustomplot.h>
#include "jahm/jahm2.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    QCustomPlot    *m_graphic;      // Declare QCustomPlot object
    QCPCurve       *m_vertical_line;     // Declare object for vertical line
    QCPCurve       *m_hyst_curve;
    ja::JAHM2       *m_jamodel;

private slots:
    void slotMousePress(QMouseEvent * event);
    void slotMouseMove(QMouseEvent * event);
};

#endif // MAINWINDOW_H
