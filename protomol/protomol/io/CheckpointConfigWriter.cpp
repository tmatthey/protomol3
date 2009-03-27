#include <protomol/io/CheckpointConfigWriter.h>

#include <iomanip>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

using namespace ProtoMol::Report;
using namespace ProtoMol;

//____CheckpointConfigWriter

CheckpointConfigWriter::CheckpointConfigWriter() :
  Writer() {}

CheckpointConfigWriter::CheckpointConfigWriter( const std::string& filename ) :
  Writer(filename) {}

bool CheckpointConfigWriter::write( const int& id, const int& steps, const Random& rand, const Integrator* integ ) {
    if (!open() ){
        return false;
    }

    file << "!Checkpoint File!" << std::endl;

    file << "#ID" << std::endl;
    file << id << std::endl;

    file << "#Step" << std::endl;
    file << steps << std::endl;

    file << "#Random" << std::endl;
    file << rand << std::endl;

    file << "#Integrator" << std::endl;
    file << (*integ) << std::endl;

    close();
    return !file.fail();
}
