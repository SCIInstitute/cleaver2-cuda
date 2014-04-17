#include "CleaverWidget.h"
#include "ui_CleaverWidget.h"
#include "MainWindow.h"
#include "ViewWidgets/MeshWindow.h"
#include <Cleaver/TetMesh.h>
#include <Cleaver/Cleaver.h>
#include <iostream>


CleaverWidget::CleaverWidget(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::CleaverWidget)
{
    ui->setupUi(this);
}

CleaverWidget::~CleaverWidget()
{
    delete ui;
}


//========================
//     Public Slots
//========================


void CleaverWidget::focus(QMdiSubWindow* subwindow)
{
    if (subwindow != NULL){

        MeshWindow *window = qobject_cast<MeshWindow *>(subwindow->widget());

        if(window != NULL)
        {
            //std::cout << "Cleaver Widget has Volume from Window" << std::endl;
            mesher = window->mesher();
        }
        else
        {
            mesher = NULL;
        }

        update();
    }
}


void CleaverWidget::clear()
{
    //
}


//=========================================
// - topologyActionChanged()
//  triggers when a topology radio button
// is selected.
//=========================================
void CleaverWidget::topologyActionChanged()
{
    if(ui->noneRadioButton->isChecked())
    {
        mesher->setTopologyMode(Cleaver::CleaverMesher::TopologyModeNone);
    }
    else if(ui->subdivideRadioButton->isChecked())
    {
        mesher->setTopologyMode(Cleaver::CleaverMesher::TopologyModeSubdivide);
    }
    else if(ui->cleaveRadioButton->isChecked())
    {
        mesher->setTopologyMode(Cleaver::CleaverMesher::TopologyModeCleave);
    }
    else
    {
        std::cout << "GUI Malfunction." << std::endl;
        exit(11);
    }
}

//=========================================
// - update()       Updates the Widget
//
//=========================================
void CleaverWidget::update()
{
    //-----------------------------
    // Set Main Mesh Button
    //-----------------------------
    if(mesher && !mesher->backgroundMeshCreated())
        ui->createMeshButton->setEnabled(true);
    else
        ui->createMeshButton->setEnabled(false);

    //-----------------------------
    // Set Background Mesh Button
    //-----------------------------
    if(mesher && mesher->backgroundMeshCreated())
        ui->createBackgroundMeshButton->setEnabled(false);
    else
        ui->createBackgroundMeshButton->setEnabled(true);

    //-----------------------------
    // Set Build Adjacency Button
    //-----------------------------
    if(mesher && mesher->backgroundMeshCreated() && !mesher->adjacencyBuilt())
        ui->buildAdjacencyButton->setEnabled(true);
    else
        ui->buildAdjacencyButton->setEnabled(false);

    //-----------------------
    //   Set Sample Button
    //-----------------------
    if(mesher && mesher->backgroundMeshCreated() && !mesher->samplingDone())
        ui->sampleVolumeButton->setEnabled(true);
    else
        ui->sampleVolumeButton->setEnabled(false);

    //----------------------------
    // Set Compute Alphas Button
    //----------------------------
    if(mesher && mesher->adjacencyBuilt() && !mesher->alphasComputed())
        ui->computeAlphasButton->setEnabled(true);
    else
        ui->computeAlphasButton->setEnabled(false);

    //-----------------------
    //    Set Cuts Button
    //-----------------------
    if(mesher && mesher->samplingDone() && !mesher->interfacesComputed())
        ui->computeInterfacesButton->setEnabled(true);
    else
        ui->computeInterfacesButton->setEnabled(false);


    //-----------------------
    // Set Generalize Button
    //-----------------------
    if(mesher && mesher->interfacesComputed() && !mesher->generalized())
        ui->generalizeTetsButton->setEnabled(true);
    else
        ui->generalizeTetsButton->setEnabled(false);

    //-----------------------
    //    Set Snap Button
    //-----------------------
    if(mesher && mesher->generalized() && !mesher->snapsAndWarpsDone())
        ui->snapAndWarpButton->setEnabled(true);
    else
        ui->snapAndWarpButton->setEnabled(false);

    //-----------------------
    //   Set Stencil Button
    //-----------------------
    if(mesher && mesher->generalized() && !mesher->stencilsDone())
        ui->stencilTetsButton->setEnabled(true);
    else
        ui->stencilTetsButton->setEnabled(false);

    /*
    //-----------------------
    //  Set Stuffing Button
    //-----------------------
    if(mesher->backgroundMesh() && !mesher->completed())
        ui->stuffButton->setEnabled(true);
    else
        ui->stuffButton->setEnabled(false);

    // allow early stenciling if generalization done
    if(mesher->generalized())
        ui->stencilButton->setEnabled(true);

    QDockWidget::update();
    */

    QDockWidget::update();
}


//=========================================
// - createMesh()
//=========================================
void CleaverWidget::createMesh()
{
    MeshWindow *window = MainWindow::instance()->activeWindow();
    if(window != NULL){
        mesher->createBackgroundMesh();

        window->setMesh(mesher->getBackgroundMesh());

        mesher->buildAdjacency();
        mesher->sampleVolume();
        window->updateMesh();
        window->updateGL();

        mesher->computeAlphas();
        mesher->computeInterfaces();
        mesher->generalizeTets();
        mesher->snapsAndWarp();

        window->updateMesh();
        window->updateGL();

        mesher->stencilTets();
        window->updateMesh();
        window->updateGL();

        Cleaver::TetMesh *mesh = mesher->getBackgroundMesh();
        mesh->name = "Adaptive-BCC-Mesh";
        MainWindow::dataManager()->addMesh(mesh);

        update();
    }
}

//=========================================
// - createBackgroundMesh()
//
//=========================================
void CleaverWidget::createBackgroundMesh()
{
    MeshWindow *window = MainWindow::instance()->activeWindow();
    if(window != NULL){
        if (!mesher->backgroundMeshCreated()){
            //Cleaver::TetMesh *mesh = Cleaver::createMeshFromVolume(window->volume());
            mesher->createBackgroundMesh();
            Cleaver::TetMesh *mesh = mesher->getBackgroundMesh();

            double total_volume = 0;
            int positive = 0;
            int negative = 0;
            for(size_t t=0; t < mesh->tets.size(); t++){
                double volume = mesh->tets[t]->volume();
                if(volume < 0)
                    negative++;
                else if(volume > 0)
                    positive++;

                total_volume += volume;
            }
            std::cout << "expected volume  = " << (double)(100*100*100) << std::endl;
            std::cout << "total tet volume = " << total_volume << std::endl;

            mesher->buildAdjacency();            

            mesh->name = "Adaptive-BCC-Mesh";
            MainWindow::dataManager()->addMesh(mesh);

            window->setMesh(mesh);
        }
        ui->computeInterfacesButton->setEnabled(true);

        //mesh->writeNodeEle("debug",true);
        update();
    }
}


//=========================================
// - buildMeshAdjacency()
//
//=========================================
void CleaverWidget::buildMeshAdjacency()
{
    mesher->buildAdjacency();
    update();
}

//=========================================
// - sampleData()
//
//=========================================
void CleaverWidget::sampleVolume()
{
    mesher->sampleVolume();
    MainWindow::instance()->activeWindow()->updateMesh();
    MainWindow::instance()->activeWindow()->updateGL();
    update();
}


//=========================================
// - computeAlphas();
//
//=========================================
void CleaverWidget::computeAlphas()
{
    mesher->computeAlphas();
    update();
}

//=========================================
// - computeCuts()
//
//=========================================
void CleaverWidget::computeInterfaces()
{
    mesher->computeInterfaces();
    MainWindow::instance()->activeWindow()->updateMesh();
    MainWindow::instance()->activeWindow()->updateGL();
    update();
}

//=========================================
// - generalizeTets()
//
//=========================================
void CleaverWidget::generalizeTets()
{
    mesher->generalizeTets();
    update();
}

//=========================================
// - snapAndWarp()
//
//=========================================
void CleaverWidget::snapAndWarp()
{
    mesher->snapsAndWarp();
    MainWindow::instance()->activeWindow()->updateMesh();
    MainWindow::instance()->activeWindow()->updateGL();
    update();
}

//=========================================
// - stencil()
//
//=========================================
void CleaverWidget::stencilTets()
{
    mesher->stencilTets();
    MainWindow::instance()->activeWindow()->updateMesh();
    MainWindow::instance()->activeWindow()->updateGL();
    update();
    mesher->getTetMesh()->writeNodeEle("output",true);
    mesher->getTetMesh()->writePly("output", true);
    //mesher->getTetMesh()->writeMatlab("mesh", true);
}

