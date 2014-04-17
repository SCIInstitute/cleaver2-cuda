#ifndef TESTDATAWIDGET_H
#define TESTDATAWIDGET_H

#include <QDockWidget>
#include <Cleaver/Volume.h>

namespace Ui {
class TestDataWidget;
}

class TestDataWidget : public QDockWidget
{
    Q_OBJECT
    
public:
    explicit TestDataWidget(QWidget *parent = 0);
    ~TestDataWidget();

public slots:

    Cleaver::Volume* createOldData();
    Cleaver::Volume* createTorusData();
    Cleaver::Volume* createConstantData();
    Cleaver::Volume* createHighDensityPatches();
    Cleaver::Volume* create2DSliceData();
    Cleaver::Volume* create2DVariableSliceData();
    Cleaver::Volume* createSmoothDensities();

    void loadTestData();
    void saveTestData();
    
private:
    Ui::TestDataWidget *ui;
};

#endif // TESTDATAWIDGET_H
