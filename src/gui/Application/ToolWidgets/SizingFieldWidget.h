#ifndef DATALOADERWIDGET_H
#define DATALOADERWIDGET_H

#include <QDockWidget>
#include <QMdiSubWindow>
#include <Cleaver/Volume.h>

namespace Ui {
class SizingFieldWidget;
}

class SizingFieldWidget : public QDockWidget
{
    Q_OBJECT
    
public:
    explicit SizingFieldWidget(QWidget *parent = 0);
    ~SizingFieldWidget();

public slots:

    void focus(QMdiSubWindow *);
    void loadIndicatorFunctions();
    void loadSizingField();
    void computeSizingField();
    void updateVolumeList();
    void volumeSelected(int index);

protected:

    void dragEnterEvent(QDragEnterEvent *event);
    void dragLeaveEvent(QDragLeaveEvent *event);
    void dropEvent(QDropEvent *event);
    
private:
    Ui::SizingFieldWidget *ui;
    Cleaver::Volume  *volume;
};

#endif // DATALOADERWIDGET_H
