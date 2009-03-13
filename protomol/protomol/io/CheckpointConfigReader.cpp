#include <protomol/io/CheckpointConfigReader.h>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

using namespace ProtoMol::Report;
using namespace ProtoMol;
//____CheckpointConfigReader

CheckpointConfigReader::CheckpointConfigReader() : Reader() {}

CheckpointConfigReader::CheckpointConfigReader(const string &filename) :
  Reader(filename) {}

bool CheckpointConfigReader::tryFormat() {
    if ( !open() ){
        return false;
    }

    std::string header;
    std::getline( file, header );
    
    std::cout << header << std::endl;

    return ( header == "!Configuration File!" );
}

bool CheckpointConfigReader::readBase( Configuration& conf, Random &rand ) {
    if (!tryFormat()) {
        std::cout << "Invalid Header" << std::endl;
        return false;
    }
    
    int id = 0, step = 0;
    
    std::string line;
    while( std::getline( file, line ) ) {
        if ( line.find("#ID") != std::string::npos ) {
            file >> id;
        }
        if ( line.find("#Step") != std::string::npos ) {
            file >> step;
        }
        if ( line.find("#Random") != std::string::npos ) {
            file >> rand;
        }
    }

    /* Update initial checkpoint perameters */
    conf["CheckpointStart"] = id;

    /* Update position file */
    conf["posfile"] = Append( conf["CheckpointPosBase"], id ) + ".pos";

    /* Update velocities file */
    conf["velfile"] = Append( conf["CheckpointVelBase"], id ) + ".vel";

    /* Update dcd file */
    if ( conf.valid( "DCDfile" ) ){
        std::string dcd = conf["DCDfile"];
        std::string dcdFile = dcd.substr( 0, dcd.rfind("dcd") + 1 );

        conf["DCDfile"] = Append( dcdFile, id ) + ".dcd";
    }

    /* Update energy file */
    if ( conf.valid( "allEnergiesFile" ) ){
        std::string energies = conf["allEnergiesFile"];

        conf["allEnergiesFile"] = Append( energies, id );
    }

    /* Update firststep */
    conf["firststep"] = toString (
        toInt( conf["firststep"] ) + step
    );

    /* Update total steps */
    conf["numsteps"] = toString (
        toInt( conf["numsteps"] ) - toInt( conf["firststep"] )
    );

    return !file.fail();
}

bool CheckpointConfigReader::readIntegrator( Integrator* integ ) {
    if (!tryFormat()) {
        return false;
    }
    
    std::string line;
    while( std::getline( file, line ) ) {
        if ( line.find("#Integrator") != std::string::npos ) {
            file >> (*integ);
        }
    }

    return !file.fail();
}
