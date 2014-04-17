// Project Includes
#include "util.h"
#include "TestData/BlobbyField.h"
#include <Cleaver/CleaverMesher.h>

// STL Includes
#include <cstdlib>

// Timing Includes
#ifdef WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

//-------------------
// Helper Function
//-------------------

float urandf()
{
    #ifdef RAND_MAX
    return ((float)rand() / RAND_MAX);
    #else
    return ((float)rand() / std::numeric_limits<int>::max());
    #endif
}

int64_t GetTime()
{
#ifdef WIN32
 /* Windows */
 FILETIME ft;
 LARGE_INTEGER li;

 /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
  * to a LARGE_INTEGER structure. */
 GetSystemTimeAsFileTime(&ft);
 li.LowPart = ft.dwLowDateTime;
 li.HighPart = ft.dwHighDateTime;

 uint64 ret = li.QuadPart;
 ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
 ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

 return ret;
#else
 /* Linux */
 struct timeval tv;

 gettimeofday(&tv, NULL);

 u_int64_t ret = tv.tv_usec;
 /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
 ret /= 1000;

 /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
 ret += (tv.tv_sec * 1000);

 return ret;
#endif
}

//---------------------------------------------------
// This is a very dumb parsing function. It requires
// ALL parameters be included, and order matters.
//---------------------------------------------------
InputParams* parseCommandline(int argc, char *argv[])
{
    int args_expected = 5;
    bool verbose = true;

    if (argc < 6){
        std::cout << "usage: "
                  << argv[0]       // argv[0]
                  << " [dimx]"     // argv[1] - dimension x
                  << " [dimy]"     // argv[2] - dimension y
                  << " [dimz]"     // argv[3] - dimension z
                  << " [mats]"     // argv[4] - # materials
                  << " [path]"     // argv[5] - output prefix path
                  << std::endl;
        exit(0);
    }

    // print to console
    if(verbose)
    {
        std::cout << "Parsing Commandline Arguments" << std::endl;
    }

    // parse dimensions
    int dim_x = atoi(argv[1]);
    int dim_y = atoi(argv[2]);
    int dim_z = atoi(argv[3]);

    // parse number of materials
    int m     = atoi(argv[4]);

    // parse prefix path
    std::string path = std::string(argv[5]);

    // create and return input param structure
    InputParams *in = new InputParams;
    in->dims = Cleaver::vec3(dim_x, dim_y, dim_z);
    in->m    = m;
    in->path = path;

    // print parsing
    if(verbose)
    {
        std::cout << "Volume Dimensions: " << in->dims.toString() << std::endl;
        std::cout << "Material Count:    " << in->m << std::endl;
        std::cout << "Path: " << in->path << std::endl;
    }

    return in;
}


//------------------------------------------------------
// Construct Input Volume
//------------------------------------------------------
Cleaver::Volume* createInputVolume(InputParams *input)
{
    // create bounds
    Cleaver::BoundingBox bounds(Cleaver::vec3::zero, input->dims);

    // create list to store all of the materials
    std::vector<Cleaver::AbstractScalarField*> fields;

    // create background material
    Cleaver::AbstractScalarField *backField = new Cleaver::ConstantField<float>(0, bounds);
    fields.push_back(backField);

    // create material indicator functions
    for(int i=0; i < input->m; i++)
    {
        // center blobby randomly between 5% and 95% from boundaries
        float rcx = (0.9*urandf() + 0.05) * 0.5*bounds.size.x + 0.25*bounds.size.x;
        float rcy = (0.9*urandf() + 0.05) * 0.5*bounds.size.y + 0.25*bounds.size.y;
        float rcz = (0.9*urandf() + 0.05) * 0.5*bounds.size.z + 0.25*bounds.size.z;

        float rr = 8 + urandf()*20;

        BlobbyField *field = new BlobbyField(Cleaver::vec3(rcx,rcy,rcz), rr, 8, bounds);
        fields.push_back(field);
    }


    // construct the volume from fields list
    Cleaver::Volume *basevol = new Cleaver::Volume(fields);
    Cleaver::Volume *volume = Cleaver::createFloatFieldVolumeFromVolume(basevol);
    volume->setName("Blobby Volume");

    return volume;
}


Cleaver::AbstractScalarField* createSizingField(Cleaver::Volume *volume, float lipschitz)
{
    return Cleaver::SizingFieldCreator::createSizingFieldFromVolume(volume, 1/lipschitz);
}


Cleaver::TetMesh* createParticleMesh(Cleaver::Volume *volume)
{
    return NULL;

    /*
    //-----------------------------
    // Sample and Optimize Particles
    //-----------------------------
    ParticleMesher particleMesher;
    particleMesher.setVolume(volume);
    particleMesher.initCPU();

    int i = 0;
    while(!particleMesher.converged() && i < 200)
    {
        particleMesher.stepCPU();
    }
    particleMesher.cleanupCPU();

    //--------------------------
    // Tetrahedralize with Tetgen
    //--------------------------
    particleMesher.tesselateCPU();

    return particleMesher.getMesh();
    */
}


//--------------------------------------------------
// Improve Mesh With Stellar
//--------------------------------------------------
void improveMeshWithStellar(Cleaver::TetMesh *mesh)
{

}

//------------------------------------------------------
// Apply Cleaving Algorithm to Mesh
//------------------------------------------------------
void cleaveMeshWithVolume(Cleaver::TetMesh *mesh, Cleaver::Volume *volume)
{
    Cleaver::CleaverMesher mesher(volume);
    mesher.setBackgroundMesh(mesh);
    mesher.buildAdjacency();
    mesher.sampleVolume();
    mesher.computeAlphas();
    mesher.computeInterfaces();
    mesher.generalizeTets();
    mesher.snapsAndWarp();
    mesher.stencilTets();
}

/*
Cleaver::Volume* loadSphereData(int n, Cleaver::BoundingBox &bounds, vector<SphereField*> &sfields)
{

//    SphereField *f1 = new SphereField(Cleaver::vec3(48.1,64,64), 12.0, bounds);
//    SphereField *f2 = new SphereField(Cleaver::vec3(80.1,64,64), 12.0, bounds);
    Cleaver::ConstantField *f3 = new Cleaver::ConstantField(0, bounds);
//    f1->setBounds(bounds);
//    f2->setBounds(bounds);
    std::vector<Cleaver::ScalarField*> fields;
//    fields.push_back(f1);
//    fields.push_back(f2);
//    fields.push_back(f3);

    for (int i=0; i<sfields.size(); i++){
        BlobbyField *bf = new BlobbyField(bounds);
    //    fields.push_back(sfields[i]);
        int num_blob = rand()%5;
        for (int j=0; j<num_blob; j++){
            Cleaver::vec3 pos;
            bf->a =  4*sfields[i]->m_r;
            bf->b = bf->a/2.0;
            bf->c = bf->b/2.0;
            pos[0] = sfields[i]->m_cx[0] + bf->a*(rand()/(double)RAND_MAX - 0.5);
            pos[1] = sfields[i]->m_cx[1] + bf->a*(rand()/(double)RAND_MAX - 0.5);
            pos[2] = sfields[i]->m_cx[2] + bf->a*(rand()/(double)RAND_MAX - 0.5);
            bf->spheres.push_back(pos);
        }
        fields.push_back(bf);
    }
    fields.push_back(f3);
    Cleaver::Volume *basevol = new Cleaver::Volume(fields);
    Cleaver::Volume *volume = Cleaver::createFloatFieldVolumeFromVolume(basevol);
    return volume;
//    volume->setName("Sphere Data");

}
*/
