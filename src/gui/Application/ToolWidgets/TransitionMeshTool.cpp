#include "TransitionMeshTool.h"
#include "ui_TransitionMeshTool.h"
#include "MainWindow.h"
#include "cleaver_utility.h"
#include <vector>
#include <array>

TransitionMeshTool::TransitionMeshTool(QWidget *parent) :
QDockWidget(parent),
ui(new Ui::TransitionMeshTool)
{
  ui->setupUi(this);

  ui->dockWidgetContents->layout()->setMargin(4);
  ui->dockWidgetContents->layout()->setSpacing(1);
}

TransitionMeshTool::~TransitionMeshTool()
{
  delete ui;
}

void TransitionMeshTool::loadVolumeData()
{
  QStringList fileNames =
      QFileDialog::getOpenFileNames(this, tr("Import Volume NRRD"),
                                    QDir::currentPath(),
                                    tr("NRRD (*.nrrd)"));
  inputs_.clear();
  for (size_t i = 0; i < fileNames.size(); i++)
    inputs_.push_back(fileNames[i].toStdString());
}

void TransitionMeshTool::createTransitionMesh()
{
  std::array<float,3> scale =
  {{ static_cast<float>(ui->xScaleSpin->value()),
      static_cast<float>(ui->yScaleSpin->value()),
      static_cast<float>(ui->zScaleSpin->value())}};
  std::array<size_t,3> res =
  {{ static_cast<size_t>(ui->xAbsSpin->value()),
      static_cast<size_t>(ui->yAbsSpin->value()),
      static_cast<size_t>(ui->zAbsSpin->value())}};
  bool absTrue = ui->absTrue->isChecked();
  bool scaleTrue = ui->scaleTrue->isChecked();
  bool useGPU = ui->useGPU->isChecked();
  bool addAir = ui->addAir->isChecked();
  CleaverUtility program;
  std::pair<std::vector<std::array<float,3>>,std::vector<std::array<size_t,4>>>
  output = program.GetVertsFacesFromNRRD(inputs_,scale.data(),res,
                                         scaleTrue,absTrue,useGPU,addAir);
  MainWindow::dataManager()->addTansitionMesh(
      output.first,output.second,{{program.w(),program.h(),program.d()}});
  MainWindow::instance()->createWindow(output.first, output.second);
}

void TransitionMeshTool::loadMesh() {
  QString fileName =
      QFileDialog::getOpenFileName(this, tr("Import Cleaver Mesh"),
                                    QDir::currentPath(),
                                    tr("CMF (*.cmf)"));
  clock_t start = clock();
  std::ifstream in(fileName.toStdString(),  std::ios::binary);
  size_t num_verts, num_faces, w, h, d;
  binary_read(in,w);
  binary_read(in,h);
  binary_read(in,d);
  std::vector<std::array<float,3>> verts;
  std::vector<std::array<size_t,4>> faces;
  binary_read(in,num_verts);
  binary_read(in,num_faces);
  for(size_t i = 0; i < num_verts; i++){
    std::array<float,3> tmp;
    binary_read(in,tmp[0]);
    binary_read(in,tmp[1]);
    binary_read(in,tmp[2]);
    verts.push_back(tmp);
  }
  for(size_t i = 0; i < num_faces; i++){
    std::array<size_t,4> tmp;
    binary_read(in,tmp[0]);
    binary_read(in,tmp[1]);
    binary_read(in,tmp[2]);
    binary_read(in,tmp[3]);
    faces.push_back(tmp);
  }
  in.close();
  double duration = ((double)clock() - (double)start) /
      (double)CLOCKS_PER_SEC;
    std::cout << "Read mesh file: " << fileName.toStdString() <<
    "\t" << duration << " sec." << std::endl;
  MainWindow::dataManager()->addTansitionMesh(verts,faces,{{w,h,d}});
  MainWindow::instance()->createWindow(verts,faces);
}

void TransitionMeshTool::outputMesh() {
  std::string name = ui->outputName->toPlainText().toStdString();
  if(name.empty()) name = "output";
  MainWindow::dataManager()->outputTansitionMesh(name);
}
