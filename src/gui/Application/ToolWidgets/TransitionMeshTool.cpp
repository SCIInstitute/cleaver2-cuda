#include "TransitionMeshTool.h"
#include "ui_TransitionMeshTool.h"
#include "MainWindow.h"
#include "cleaver_utility.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>

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
