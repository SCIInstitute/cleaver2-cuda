#include "TestDataWidget.h"
#include "ui_TestDataWidget.h"
#include "MainWindow.h"
#include <Cleaver/Cleaver.h>
#include <Cleaver/InverseField.h>
#include <Cleaver/ConstantField.h>
#include "TestData/PlaneField.h"
#include "TestData/SphereField.h"
#include "TestData/SphereSizingField.h"
#include "TestData/SphereVaryingField.h"
#include "TestData/PlaneSizingField.h"
#include "TestData/TorusField.h"
#include "TestData/BlobbyField.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>

TestDataWidget::TestDataWidget(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::TestDataWidget)
{
    ui->setupUi(this);

    ui->dockWidgetContents->layout()->setMargin(4);
    ui->dockWidgetContents->layout()->setSpacing(1);
    ui->dockWidgetContents->layout()->setAlignment(ui->loadTestDataButton, Qt::AlignHCenter);
    ui->dockWidgetContents->layout()->setAlignment(ui->saveTestDataButton, Qt::AlignHCenter);
}

TestDataWidget::~TestDataWidget()
{
    delete ui;
}

namespace TestData{
float randf()
{
    #ifdef RAND_MAX
    return ((float)rand() / RAND_MAX);
    #else
    return ((float)rand() / std::numeric_limits<int>::max());
    #endif
}
}


Cleaver::Volume* TestDataWidget::createTorusData()
{
    std::string dataTitle("Torus Data");
    std::cout << "Loading " << dataTitle << std::endl;

    Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(64,64,64));

    std::vector<Cleaver::AbstractScalarField*> fields;

    //Cleaver::ConstantField *f0 = new Cleaver::ConstantField(0, bounds);
    //f0->setName("Constant Field");
    TorusField *f1 = new TorusField(Cleaver::vec3(32, 32, 32), 12, 5, bounds);   // breaks feature size computation
    f1->setName("Torus");
    Cleaver::AbstractScalarField *f2 = new Cleaver::InverseScalarField(f1);
    f2->setName("InverseTorus");


    fields.push_back(f1);
    fields.push_back(f2);

    Cleaver::Volume *basevol = new Cleaver::Volume(fields);
    Cleaver::Volume *volume = Cleaver::createFloatFieldVolumeFromVolume(basevol);
    volume->setName(dataTitle.c_str());

    return volume;

    //MainWindow::instance()->createWindow(basevol, dataTitle.c_str());
}

Cleaver::Volume* TestDataWidget::createConstantData()
{
    std::string dataTitle("Constant Field Test");
    std::cout << "Loading " << dataTitle << std::endl;

    int w = 32;
    int h = 32;
    int d = 32;
    Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(w,h,d));

    std::vector<Cleaver::AbstractScalarField*> fields;

    Cleaver::ConstantField<float> *f0 = new Cleaver::ConstantField<float>(0, bounds);
    TorusField *f1 = new TorusField(Cleaver::vec3(16, 16, 16), 6, 2.5, bounds);   // breaks feature size computation
    //SphereField *f2 = new SphereField(Cleaver::vec3(32,32,32), 4, bounds);

    fields.push_back(f0);
    fields.push_back(f1);

    Cleaver::Volume *basevol = new Cleaver::Volume(fields);
    basevol->setName(dataTitle.c_str());

    Cleaver::Volume *volume = Cleaver::createFloatFieldVolumeFromVolume(basevol);
    volume->setName(dataTitle.c_str());

    Cleaver::ConstantField<float> *sizingField = new Cleaver::ConstantField<float>(5, bounds);
    Cleaver::FloatField *sizingFloatField = Cleaver::createFloatFieldFromScalarField(sizingField);


    basevol->setSizingField(sizingFloatField);
    std::cout << " = " << volume->name() << std::endl;

    return basevol;
}

Cleaver::Volume* TestDataWidget::createHighDensityPatches()
{
    std::string dataTitle("Constant Field Test");
    std::cout << "Loading " << dataTitle << std::endl;

    int w = 32;
    int h = 32;
    int d = 32;
    Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(w,h,d));

    std::vector<Cleaver::AbstractScalarField*> fields;

    Cleaver::ConstantField<float> *f0 = new Cleaver::ConstantField<float>(0, bounds);
    TorusField *f1 = new TorusField(Cleaver::vec3(16, 16, 16), 6, 2.5, bounds);   // breaks feature size computation
    //SphereField *f2 = new SphereField(Cleaver::vec3(32,32,32), 4, bounds);

    fields.push_back(f0);
    fields.push_back(f1);

    Cleaver::Volume *basevol = new Cleaver::Volume(fields);
    basevol->setName(dataTitle.c_str());

    Cleaver::Volume *volume = Cleaver::createFloatFieldVolumeFromVolume(basevol);
    volume->setName(dataTitle.c_str());

    Cleaver::ConstantField<float> *sizingField = new Cleaver::ConstantField<float>(5, bounds);
    Cleaver::FloatField *sizingFloatField = Cleaver::createFloatFieldFromScalarField(sizingField);


    basevol->setSizingField(sizingFloatField);
    std::cout << " = " << volume->name() << std::endl;

    return basevol;
}


Cleaver::Volume* TestDataWidget::create2DSliceData()
{
    int w = 31;
    int h = 31;
    int d = 31;
    std::string dataTitle("2D-Constant");
    std::cout << "Loading " << dataTitle << std::endl;

    Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(w,h,d));

    std::vector<Cleaver::AbstractScalarField*> fields;

    Cleaver::ConstantField<float> *f0 = new Cleaver::ConstantField<float>(0, bounds);
    f0->setName(std::string("background"));
    TorusField *f1 = new TorusField(Cleaver::vec3(w/2, h/2, d/2), 6, 2.5, bounds);   // breaks feature size computation
    f1->setName("torus");
    //SphereField *f2 = new SphereField(Cleaver::vec3(32,32,32), 4, bounds);

    fields.push_back(f0);
    fields.push_back(f1);

    Cleaver::Volume *basevol = new Cleaver::Volume(fields);
    Cleaver::Volume *volume = Cleaver::createFloatFieldVolumeFromVolume(basevol);

    basevol->setName(dataTitle.c_str());
    volume->setName(dataTitle.c_str());

    float *data = new float[w*h*d];

    // set background to be effectively empty
    for(int k=0; k < d; k++)
    {
        for(int j=0; j < h; j++)
        {
            for(int i=0; i < w; i++)
            {
                data[k*(w*h) + j*w + i] = 1e8;
            }
        }
    }

    // set sizing field on the one central plane we care about
    int k = d / 2;
    for(int j=0; j < h; j++)
    {
        for(int i=0; i < w; i++)
        {
            data[k*(w*h) + j*w + i] = 2;
        }
    }

    Cleaver::FloatField *sizingFloatField = new Cleaver::FloatField(data, w, h, d);
    sizingFloatField->setName("constant-plane");
    MainWindow::dataManager()->addField(sizingFloatField);


    basevol->setSizingField(sizingFloatField);
    return basevol;
}

Cleaver::Volume* TestDataWidget::create2DVariableSliceData()
{
    int w = 31;
    int h = 31;
    int d = 31;
    std::string dataTitle("2D Constant Field Test");
    std::cout << "Loading " << dataTitle << std::endl;

    Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(w,h,d));

    std::vector<Cleaver::AbstractScalarField*> fields;

    Cleaver::ConstantField<float> *f0 = new Cleaver::ConstantField<float>(0, bounds);
    f0->setName(std::string("background"));
    TorusField *f1 = new TorusField(Cleaver::vec3(w/2, h/2, d/2), 6, 2.5, bounds);   // breaks feature size computation
    f1->setName(std::string("torus"));


    fields.push_back(f0);
    fields.push_back(f1);

    Cleaver::Volume *basevol = new Cleaver::Volume(fields);
    Cleaver::Volume *volume = Cleaver::createFloatFieldVolumeFromVolume(basevol);
    basevol->setName(dataTitle.c_str());
    volume->setName(dataTitle.c_str());

    float *data = new float[w*h*d];

    // set background to be effectively empty
    for(int k=0; k < d; k++)
    {
        for(int j=0; j < h; j++)
        {
            for(int i=0; i < w; i++)
            {
                data[k*(w*h) + j*w + i] = 1e8;
            }
        }
    }

    // set sizing field on the one central plane we care about
    int k = d / 2;
    for(int j=0; j < h; j++)
    {
        for(int i=0; i < w; i++)
        {
            float cx = w/2.0; float cy = h/2.0;
            float r = sqrt((cx-i)*(cx-i) + (cy-j)*(cy-j));
            data[k*(w*h) + j*w + i] = 1.1 + 3*(i/31.0);
            //std::cout << r << std::endl;
        }
    }

//    data[15*(w*h) + 15*w + 15] = 0.5;

    Cleaver::FloatField *sizingFloatField = new Cleaver::FloatField(data, w, h, d);


    basevol->setSizingField(sizingFloatField);
    return basevol;
}

Cleaver::Volume* TestDataWidget::createOldData()
{
    std::cout << "Loading Old Test Data" << std::endl;

    Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(64,64,64));
    //Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(32,32,32));

    //Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(36,36,36));
    //Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(16,16,16));


    //SphereField *f1 = new SphereField(Cleaver::vec3(26.5,32,32), 7.2, bounds);
    //SphereField *f2 = new SphereField(Cleaver::vec3(38.5,32,32), 7.2, bounds);

    //PlaneField *f1 = new PlaneField(Cleaver::vec3(0.3,+1.0,0.2), Cleaver::vec3(32,32.4,32));   // VERY BAD #1
    //PlaneField *f2 = new PlaneField(Cleaver::vec3(0.4,-1.0,0), Cleaver::vec3(32,32.4,32));     // VERY BAD #1

    //PlaneField *f1 = new PlaneField(Cleaver::vec3(0.3,+1.0,0.2), Cleaver::vec3(32,32.4,32));     // worst +0.5,+1.0,+0.5
    //PlaneField *f2 = new PlaneField(Cleaver::vec3(0.4,-1.0,0), Cleaver::vec3(32,32.4,32));     // worst +0.5,-1.0,+0.5

    //PlaneField *f1 = new PlaneField(Cleaver::vec3(0.3,+1.0,0.2), Cleaver::vec3(16,16.2,16));     // worst +0.5,+1.0,+0.5
    //PlaneField *f2 = new PlaneField(Cleaver::vec3(0.4,-1.0,0), Cleaver::vec3(16,16.2,16));     // worst +0.5,-1.0,+0.5

    //SphereField *f1 = new SphereField(Cleaver::vec3(8.1, 8 , 8), 4.0, bounds);
    //SphereField *f2 = new SphereField(Cleaver::vec3(20.1, 22, 20), 8.0, bounds);
//    SphereField *f1 = new SphereField(Cleaver::vec3(8.1+16, 8+16 , 8+16), 4.0, bounds);
//    SphereField *f2 = new SphereField(Cleaver::vec3(20.1+16, 22+16, 20+16), 8.0, bounds);

//    SphereField *f1 = new SphereField(Cleaver::vec3(24.1, 24 , 24), 4.0, bounds);
//    SphereField *f2 = new SphereField(Cleaver::vec3(20.1+16, 22+16, 20+16), 8.0, bounds);

    //TorusField *f1 = new TorusField(Cleaver::vec3(17.5f,17.5f,17.5f), 8, 3, bounds);   // breaks feature size computation


    //Cleaver::InverseField *f2 = new Cleaver::InverseField(f1);

    //SphereField *f1 = new SphereField(Cleaver::vec3(16,16,16), 10, bounds);

    SphereField *f3 = new SphereField(Cleaver::vec3(32,32,32), 4, bounds);

    std::vector<Cleaver::AbstractScalarField*> fields;
    Cleaver::ConstantField<float> *f0 = new Cleaver::ConstantField<float>(0.0f, bounds);
    //SphereField *f2 = new SphereField(Cleaver::vec3(32,32,32), 20, bounds);
    TorusField *f1 = new TorusField(Cleaver::vec3(32, 32, 32), 12, 5, bounds);   // breaks feature size computation
    //Cleaver::InverseField *f2 = new Cleaver::InverseField(f1);

    fields.push_back(f0);
    fields.push_back(f1);

    //fields.push_back(f2);
    //fields.push_back(f2);
    //fields.push_back(f3);

    Cleaver::Volume *basevol = new Cleaver::Volume(fields);
    Cleaver::Volume *volume = Cleaver::createFloatFieldVolumeFromVolume(basevol);
    basevol->setName("Sphere Data");
    volume->setName("Sphere Data");

    //delete f1;
    //delete f2;
    //delete f3;

    //SphereSizingField *sizingField = new SphereSizingField(f1);
    //SphereVaryingField *sizingField2 = new SphereVaryingField(f1);

    //volume->setSizingField(sizingField);
    //volume->setSizingField(sizingField2);


    float *fdata = new float[(int)(bounds.size.x * bounds.size.y * bounds.size.z)];
    memset(fdata, 0, bounds.size.x * bounds.size.y * bounds.size.z * sizeof(float));
    int idx=0;
    for(int k=0; k < bounds.size.z; k++)
    {
        for(int j=0; j < bounds.size.y; j++)
        {
            for(int i=0; i < bounds.size.x; i++)
            {
                fdata[idx++] = 16;
            }
        }
    }
    int pos[3]= {6,6,6};


    int w = bounds.size.x;
    int h = bounds.size.y;
    int d = bounds.size.z;

    for(int k=10; k < 13; k++)
    {
        for(int j=10; j < 13; j++)
        {
            for(int i=10; i < 13; i++)
            {
                fdata[k*w*h + j*w + i] = 0.4;
            }
        }
    }

    Cleaver::FloatField *fieldTest = new Cleaver::FloatField(fdata, bounds.size.x, bounds.size.y, bounds.size.z);
    basevol->setSizingField(fieldTest);

    //Cleaver::ConstantField *constantField = new Cleaver::ConstantField(5, bounds);
    //volume->setSizingField(constantField);

    //Cleaver::ConstantField *constantField = new Cleaver::ConstantField(8, Cleaver::BoundingBox(Cleaver::vec3(0,0,0), Cleaver::vec3(64,64,64)));
    //Cleaver::FloatField *constantFloatField = Cleaver::createFloatFieldFromScalarField(constantField);
    //constantFloatField->setBounds(Cleaver::BoundingBox(Cleaver::vec3(-16,-16,-16), Cleaver::vec3(64,64,64)));
    //constantFloatField->setBounds(Cleaver::BoundingBox(Cleaver::vec3(-16,-16,-16), Cleaver::vec3(64,64,64)));
    //volume->setSizingField(constantFloatField);

//    QString title = "Sphere Data";
//    MainWindow::instance()->createWindow(basevol, title);
    return basevol;
}

Cleaver::Volume* TestDataWidget::createSmoothDensities()
{
    Cleaver::BoundingBox bounds(Cleaver::vec3::zero, Cleaver::vec3(64,64,64));

    int w = bounds.size.x;
    int h = bounds.size.y;
    int d = bounds.size.z;

    Cleaver::ConstantField<float> *f0 = new Cleaver::ConstantField<float>(0.0f, bounds);
    TorusField *f1 = new TorusField(Cleaver::vec3(32, 32, 32), 12, 5, bounds);   // breaks feature size computation

    std::vector<Cleaver::AbstractScalarField*> fields;
    fields.push_back(f0);
    fields.push_back(f1);

    float *fdata = new float[(int)(bounds.size.x * bounds.size.y * bounds.size.z)];
    memset(fdata, 0, bounds.size.x * bounds.size.y * bounds.size.z * sizeof(float));


    Cleaver::vec3 c1(8, 8, 8);
    Cleaver::vec3 c2(32, 53, 42);

    int idx=0;
    for(int k=0; k < bounds.size.z; k++)
    {
        for(int j=0; j < bounds.size.y; j++)
        {
            for(int i=0; i < bounds.size.x; i++)
            {
                float r1 = 0.5*std::abs(Cleaver::L2(Cleaver::vec3(i,j,k) - c1));
                float r2 = 0.5*std::abs(Cleaver::L2(Cleaver::vec3(i,j,k) - c2));

                fdata[idx++] = std::max(std::min(std::min(16.0f, r1), r2), 2.0f);
            }
        }
    }

    Cleaver::FloatField *fieldTest = new Cleaver::FloatField(fdata, bounds.size.x, bounds.size.y, bounds.size.z);

    Cleaver::Volume *volume = new Cleaver::Volume(fields);
    volume->setName("Smooth Densities");
    volume->setSizingField(fieldTest);

    return volume;
}

void TestDataWidget::loadTestData()
{
    Cleaver::Volume *volume = createTorusData();
    //Cleaver::Volume *volume = create2DSliceData();
    //Cleaver::Volume *volume = create2DVariableSliceData();
    //Cleaver::Volume *volume = createConstantData();
    //Cleaver::Volume *volume = createHighDensityPatches();
    //Cleaver::Volume *volume = createOldData();
    //Cleaver::Volume *volume = createSmoothDensities();

    for(size_t f=0; f < volume->numberOfMaterials(); f++)
        MainWindow::dataManager()->addField(volume->getMaterial(f));

    MainWindow::dataManager()->addVolume(volume);
    MainWindow::instance()->createWindow(volume, volume->name().c_str());
}

void TestDataWidget::saveTestData()
{
    MeshWindow *window = MainWindow::instance()->activeWindow();
    if(!window)
        return;

    Cleaver::Volume *volume = window->volume();
    if(!volume)
        return;

    for(int i=0; i < volume->numberOfMaterials(); i++)
    {
        Cleaver::FloatField *field = (Cleaver::FloatField*)volume->getMaterial(i);

        std::stringstream ss;

        ss << "material_" << i << ".nrrd";
        std::string filename = ss.str();

        int w = field->bounds().size.x;
        int h = field->bounds().size.y;
        int d = field->bounds().size.z;

        std::cout << "Writing file " << filename << std::endl;
        std::ofstream nrrd_file(filename.c_str());
        nrrd_file << "NRRD0004" << std::endl;
        nrrd_file << "# Complete NRRD file format specification at:" << std::endl;
        nrrd_file << "# http://teem.sourceforge.net/nrrd/format.html" << std::endl;
        nrrd_file << "type: float" << std::endl;
        nrrd_file << "dimension: 3" << std::endl;
        nrrd_file << "space dimension: 3" << std::endl;
        nrrd_file << "sizes: " << w << " " << h << " " << d << std::endl;
        //nrrd_file << "space directions: (1,0,0) (0,1,0) (0,0,1)" << endl;
        //nrrd_file << "space directions: (" << sx << ",0,0) (0," << sy << ",0) (0,0," << sz << ")" << std::endl;
        nrrd_file << "centerings: node node node" << std::endl;
        nrrd_file << "endian: little" << std::endl;
        nrrd_file << "encoding: raw" << std::endl;
        nrrd_file << "space origin: (0,0,0)" << std::endl << std::endl;

            // write data portion
            for(int k=0; k < d; k++)
            {
                for(int j=0; j < h; j++)
                {
                    for(int i=0; i < w; i++)
                    {
                        float val = field->valueAt(i+0.5, j+0.5, k+0.5);
                        nrrd_file.write((char*)&val, sizeof(float));
                    }
                }
            }

            nrrd_file.close();
    }
}
