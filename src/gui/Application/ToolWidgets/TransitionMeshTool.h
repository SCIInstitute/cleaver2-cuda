#ifndef TRANSITIONMESHTOOL_H
#define TRANSITIONMESHTOOL_H

#include <QDockWidget>

namespace Ui {
class TransitionMeshTool;
}

class TransitionMeshTool : public QDockWidget
{
    Q_OBJECT
    
public:
    explicit TransitionMeshTool(QWidget *parent = 0);
    ~TransitionMeshTool();

public slots:

    void loadVolumeData();
    void createTransitionMesh();
    
private:
    Ui::TransitionMeshTool *ui;
    std::vector<std::string> inputs_;
};

#endif // TRANSITIONMESHTOOL_H
