// Utils Includes
#include "util.h"

// STL Includes
#include <iostream>
#include <ctime>
#include <sstream>
#include <fstream>
#include <string>


// Entry Point
int main(int argc,	char* argv[])
{
    int64_t begin_time = GetTime();

    //-----------------------------------------------------------
    // Parse Command-line Params
    //-----------------------------------------------------------
    InputParams *input = parseCommandline(argc, argv);


    //-----------------------------------------------------------
    // Construct Input Volume
    //-----------------------------------------------------------
    Cleaver::Volume *volume = createInputVolume(input);


    //------------------------------------------------------------
    // Construct Sizing Field
    //------------------------------------------------------------
    Cleaver::AbstractScalarField *sizingField = createSizingField(volume, 0.2);
    volume->setSizingField(sizingField);


    //-----------------------------------------------------------
    // Construct Particle Mesh
    //-----------------------------------------------------------
    Cleaver::TetMesh *mesh = createParticleMesh(volume);

    // write to file
    std::string particleMeshName = input->path + "/particlemesh";
    mesh->writeNodeEle(particleMeshName, false, false, false);
    mesh->writePly(particleMeshName);
    mesh->writeInfo(particleMeshName);

    //-----------------------------------------------------------
    // Improve With Stellar
    //-----------------------------------------------------------
    improveMeshWithStellar(mesh);

    // write to file
    std::string stellarMeshName = input->path + "/stellarmesh";
    mesh->writeNodeEle(stellarMeshName, false, false, false);
    mesh->writePly(stellarMeshName);
    mesh->writeInfo(stellarMeshName);


    //-----------------------------------------------------------
    // Apply Mesh Cleaving
    //-----------------------------------------------------------
    cleaveMeshWithVolume(mesh, volume);

    // write to file
    std::string cleaverMeshName = input->path + "/cleavermesh";
    mesh->writeNodeEle(cleaverMeshName, false, true, true);
    mesh->writePly(cleaverMeshName);
    mesh->writeInfo(cleaverMeshName);

    //-----------------------------------------------------------
    // Write Experiment Info to file
    //-----------------------------------------------------------
    int64_t end_time = GetTime();
    double total_time = (end_time - begin_time)/(double)1000;


    std::string infoFilename = input->path + "/experiment.info";
    std::ofstream file(infoFilename.c_str());
    file << "Experiment Info" << std::endl;
    file << "Size: " << input->dims.toString() << std::endl;
    file << "Materials: " << input->m << std::endl;
    file << "Time: " << total_time << " seconds" << std::endl;

}
