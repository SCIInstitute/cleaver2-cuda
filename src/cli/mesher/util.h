//------------------------------------------------------
// Cleaver Includes
#include <Cleaver/CleaverMesher.h>
#include <Cleaver/ScalarField.h>
#include <Cleaver/ConstantField.h>
#include <Cleaver/Cleaver.h>
#include <Cleaver/SizingFieldCreator.h>

#include <string>

//------------------------------------------------------
// Command-line Parameters
//------------------------------------------------------
struct InputParams
{
    Cleaver::vec3 dims;        // volume dimensions
    int m;                     // number of materials
    std::string path;    // output mesh filename
};



//------------------------------------------------------
// Parse Command-line Params
//------------------------------------------------------
InputParams* parseCommandline(int argc, char* argv[]);


//------------------------------------------------------
// Construct Input Volume
//------------------------------------------------------
Cleaver::Volume* createInputVolume(InputParams *input);


//------------------------------------------------------
// Construct Sizing Field
//------------------------------------------------------
Cleaver::AbstractScalarField* createSizingField(Cleaver::Volume *volume, float lipschitz);

//------------------------------------------------------
// Construct Particle Mesh
//------------------------------------------------------
Cleaver::TetMesh* createParticleMesh(Cleaver::Volume *volume);

//------------------------------------------------------
// Improve Mesh With Stellar
//------------------------------------------------------
void improveMeshWithStellar(Cleaver::TetMesh *mesh);

//------------------------------------------------------
// Apply Cleaving Algorithm to Mesh
//------------------------------------------------------
void cleaveMeshWithVolume(Cleaver::TetMesh *mesh, Cleaver::Volume *volume);


int64_t GetTime();

//Cleaver::Volume* loadSphereData(int n, Cleaver::BoundingBox &bounds, std::vector<Cleaver::SphereField*> &sfields);
