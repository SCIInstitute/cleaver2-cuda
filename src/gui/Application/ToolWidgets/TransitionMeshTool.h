#ifndef TRANSITIONMESHTOOL_H
#define TRANSITIONMESHTOOL_H

#include <QDockWidget>
#include <iostream>
#include <fstream>
#include <Cleaver/Cleaver.h>

template<typename T>
std::istream & binary_read(std::istream& stream, T& value){
    return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}
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
    void loadMesh();
    void outputMesh();
    
private:
    Ui::TransitionMeshTool *ui;
    std::vector<std::string> inputs_;
    Cleaver::FloatField * ff_;
};

#endif // TRANSITIONMESHTOOL_H
