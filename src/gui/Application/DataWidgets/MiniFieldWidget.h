#ifndef MINIFIELDWIDGET_H
#define MINIFIELDWIDGET_H

#include <QWidget>
#include <Cleaver/ScalarField.h>

namespace Ui {
class MiniFieldWidget;
}

class MiniFieldWidget : public QWidget
{
    Q_OBJECT
    
public:
    explicit MiniFieldWidget(QWidget *parent = 0);
    MiniFieldWidget(Cleaver::AbstractScalarField *field, QWidget *parent = 0);
    ~MiniFieldWidget();

    void setField(Cleaver::AbstractScalarField *field);
    Cleaver::AbstractScalarField* field() const {  return m_field; }

signals:

    void removeRequest(Cleaver::AbstractScalarField*);
    void changeRequest(Cleaver::AbstractScalarField*);

protected:

    void paintEvent(QPaintEvent *);

    void mousePressEvent(QMouseEvent *);
    void mouseReleaseEvent(QMouseEvent *);
    void mouseMoveEvent(QMouseEvent *);


    void startDrag(Qt::DropActions supportedActions);
    void dragEnterEvent(QDragEnterEvent *event);
    void dragLeaveEvent(QDragLeaveEvent *event);
    void dropEvent(QDropEvent *event);

    
private:
    Ui::MiniFieldWidget *ui;
    Cleaver::AbstractScalarField *m_field;

    std::string normalStyle;
    std::string dragEnterStyle;
};

#endif // MINIFIELDWIDGET_H
