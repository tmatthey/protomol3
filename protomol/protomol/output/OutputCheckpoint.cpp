#include <protomol/output/OutputCheckpoint.h>

#include <protomol/topology/TopologyUtilities.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/base/MathUtilities.h>

#include <protomol/io/XYZWriter.h>

#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

const string OutputCheckpoint::keyword( "Checkpoint" );

OutputCheckpoint::OutputCheckpoint() : mName( "false" ), mCurrent( 0 ) {
    if ( mName == "true" || mName == "True" || mName == "TRUE" ) {
        isActive = true;
    }else{
        isActive = false;
    }
}

OutputCheckpoint::OutputCheckpoint( const std::string& name, int freq, int start ) :
        Output( freq ), mName( name ), mCurrent( start ) {
    if ( mName == "true" || mName == "True" || mName == "TRUE" ) {
        isActive = true;
    }else{
        isActive = false;
    }
}

void OutputCheckpoint::doInitialize() {
    if ( isActive ) {
        const std::string temp = app->config["posfile"];

        mBaseFile = temp.substr( 0, temp.rfind( '.' ) );

        std::cout << "Checkpointing: Active" << std::endl;
        std::cout << "Checkpointing: Base " << mBaseFile << std::endl;
    }
}

void OutputCheckpoint::doRun( int step ) {
    if ( step != 0 ) {
        if ( isActive ) {
            ReadConfig();
            WritePositions( step );
            WriteVelocities( step );
            WriteConfig( step );

            mCurrent += 1;

            std::cout << "Checkpointing: Step " << step << std::endl;
        }
    }
}

void OutputCheckpoint::doFinalize( int step ) {

}

Output *OutputCheckpoint::doMake( const vector<Value> &values ) const {
    return new OutputCheckpoint( values[0], values[1], values[2] );
}

bool OutputCheckpoint::isIdDefined( const Configuration *config ) const {
    return config->valid( getId() );
}

void OutputCheckpoint::getParameters( vector<Parameter> &parameter ) const {
    parameter.push_back(
        Parameter( getId(),
                   Value( mName, ConstraintValueType::NotEmpty() ) ) );

    parameter.push_back(
        Parameter( keyword + "Freq",
                   Value( myOutputFreq, ConstraintValueType::Positive() ) ) );

    parameter.push_back(
        Parameter( keyword + "Start",
                   Value( mCurrent, ConstraintValueType::NotNegative() ) ) );
}

bool OutputCheckpoint::adjustWithDefaultParameters( vector<Value> &values,
        const Configuration *config ) const {
    if ( !checkParameterTypes( values ) ) return false;

    if ( !values[0].valid() ) {
        values[0] = mName;
    }

    if ( !values[1].valid() ) {
        values[1] = 1;
    }

    if ( !values[2].valid() ) {
        values[2] = 0;
    }

    return checkParameters( values );
}

void OutputCheckpoint::ReadConfig( ) {
    ifstream inFile( "checkpoint.dat" );
    ofstream outFile( "checkpoint.last" );

    if ( inFile && outFile ) {
        std::string line;

        while ( std::getline( inFile, line ) ) {
            outFile << line << '\n';
        }
    }
}

void OutputCheckpoint::WritePositions( int step ) {
    std::string posFile = Append( mBaseFile + ".pos", mCurrent );

    XYZWriter posWriter;
    if ( !posWriter.open( posFile ) ) {
        THROW( string( "Can't open " ) + getId() + " '" + posFile + "'." );
    }

    const Vector3DBlock *pos = &app->positions;
    posWriter.setComment( "Time : " + toString( app->outputCache.time() ) +
                          ", step : " + toString( step ) +  "." );

    if ( !posWriter.write( *pos, app->topology->atoms, app->topology->atomTypes ) ) {
        THROW( string( "Could not write " ) + getId() + " '" + posFile + "'." );
    }
}

void OutputCheckpoint::WriteVelocities( int step ) {
    std::string velFile = Append( mBaseFile + ".vel", mCurrent );

    XYZWriter velWriter;
    if ( !velWriter.open( velFile ) ) {
        THROW( string( "Can't open " ) + getId() + " '" + velFile + "'." );
    }

    velWriter.setComment( "Time : " + toString( app->outputCache.time() ) +
                          ", step : " + toString( step ) + "." );

    if ( !velWriter.write( *&app->velocities, app->topology->atoms,
                           app->topology->atomTypes ) ) {
        THROW( string( "Could not write " ) + getId() + " '" + velFile + "'." );
    }
}

void OutputCheckpoint::WriteConfig( int step ){
    ofstream outFile( "checkpoint.dat" );

    if ( outFile ) {
        outFile << mCurrent << std::endl;
        outFile << step << std::endl;
        outFile << Random::Instance() << std::endl;
    }
}
