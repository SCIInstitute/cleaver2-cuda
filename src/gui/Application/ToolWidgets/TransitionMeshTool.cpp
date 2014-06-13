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
      QFileDialog::getOpenFileName(this, tr("Import Polygon Mesh"),
                                   QDir::currentPath(),
                                   tr("PLY (*.ply)"));
  clock_t start = clock();
  std::ifstream in(fileName.toStdString());
  std::string tmp;
  size_t num_verts=0, num_faces=0, width=0, height=0, depth=0;
  unsigned long n_verts=0, n_faces=0, w=0, h=0, d=0;
  getline(in,tmp);
  getline(in,tmp);
  getline(in,tmp);
  sscanf(tmp.c_str(),"comment %lu %lu %lu",
         &w, &h, &d);
  width = static_cast<size_t>(w);
  height = static_cast<size_t>(h);
  depth = static_cast<size_t>(d);
  getline(in,tmp);
  sscanf(tmp.c_str(),"element vertex %lu",&n_verts);
  getline(in,tmp);
  getline(in,tmp);
  getline(in,tmp);
  getline(in,tmp);
  sscanf(tmp.c_str(),"element face %lu",&n_faces);
  getline(in,tmp);
  getline(in,tmp);
  getline(in,tmp);
  getline(in,tmp);
  num_verts = static_cast<size_t>(n_verts);
  num_faces = static_cast<size_t>(n_faces);
  std::vector<std::array<float,3>> verts;
  std::vector<std::array<size_t,4>> faces;
  for(size_t i = 0; i < num_verts; i++){
    getline(in,tmp);
    std::array<float,3> v;
    sscanf(tmp.c_str(),"%g %g %g",&v[0],&v[1],&v[2]);
    verts.push_back(v);
  }
  size_t dum;
  for(size_t i = 0; i < num_faces; i++){
    getline(in,tmp);
    std::array<size_t,4> f;
    std::array<unsigned long,4> ff;
    sscanf(tmp.c_str(),"%lu %lu %lu %lu",&dum,&ff[0],&ff[1],&ff[2]);

    for(size_t j = 0; j < 3; j++) f[j] = static_cast<size_t>(ff[j]);
    faces.push_back(f);
  }
  for(size_t i = 0; i < num_faces; i++){
    getline(in,tmp);
    unsigned long nn = 0;
    sscanf(tmp.c_str(),"%lu",&nn);
    faces.at(i)[3] = static_cast<size_t>(nn);
  }
  in.close();
  double duration = ((double)clock() - (double)start) /
      (double)CLOCKS_PER_SEC;
  std::cout << "Read mesh file: " << fileName.toStdString() <<
      "\t" << duration << " sec." << std::endl;
  MainWindow::dataManager()->addTansitionMesh(
      verts,faces,{{width,height,depth}});
  MainWindow::instance()->createWindow(verts,faces);
}

void TransitionMeshTool::outputMesh() {
  std::string name = ui->outputName->toPlainText().toStdString();
  if(name.empty()) name = "output";
  MainWindow::dataManager()->outputTansitionMesh(name);
}
