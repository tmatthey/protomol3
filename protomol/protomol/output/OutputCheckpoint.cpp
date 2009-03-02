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

OutputCheckpoint::OutputCheckpoint() : mName( "false" ), mCurrent( 0 ),
        mPosBase( "" ), mVelBase( "" ) {
    if ( mName == "true" || mName == "True" || mName == "TRUE" ) {
        isActive = true;
    } else {
        isActive = false;
    }
}

OutputCheckpoint::OutputCheckpoint( const std::string& name, int freq,
                                    int start, const std::string& posbase,
                                    const std::string& velbase ) :
        Output( freq ), mName( name ), mCurrent( start ),
        mPosBase( posbase ), mVelBase( velbase ) {

    if ( mName == "true" || mName == "True" || mName == "TRUE" ) {
        isActive = true;
    } else {
        isActive = false;
    }
}

void OutputCheckpoint::doInitialize() {
    if ( isActive ) {
        if ( mPosBase == "" ) {
            const std::string temp = app->config["posfile"];

            mPosBase = temp.substr( 0, temp.rfind( '.' ) + 1 );
        }

        if ( mVelBase == "" ) {
            if ( !app->config.valid( "velfile" ) ){
                mVelBase = mPosBase;
            }else{
                const std::string temp = app->config["velfile"];

                mVelBase = temp.substr( 0, temp.rfind( '.' ) + 1 );
            }
        }
    }
}

void OutputCheckpoint::doRun( int step ) {
    if ( isActive ) {
        if ( step != 0 ){
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
    return new OutputCheckpoint( values[0], values[1], values[2], values[3],
                                 values[4] );
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

    parameter.push_back(
        Parameter( keyword + "PosBase",
                   Value( mPosBase, ConstraintValueType::NoConstraints() ) ) );

    parameter.push_back(
        Parameter( keyword + "VelBase",
                   Value( mVelBase, ConstraintValueType::NoConstraints() ) ) );
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

    if ( !values[3].valid() ) {
        values[3] = "";
    }

    if ( !values[4].valid() ) {
        values[4] = "";
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
    std::string posFile = Append( Append( mPosBase, mCurrent ), ".pos" );

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
    std::string velFile = Append( Append( mVelBase, mCurrent ), ".vel" );

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

void OutputCheckpoint::WriteConfig( int step ) {
    ofstream outFile( "checkpoint.dat" );

    if ( outFile ) {
        outFile << "#ID" << std::endl;
        outFile << mCurrent << std::endl;

        outFile << "#Step" << std::endl;
        outFile << step << std::endl;

        outFile << "#Random" << std::endl;
        outFile << Random::Instance() << std::endl;

        outFile << "#Integrator" << std::endl;
        outFile << (*app->integrator) << std::endl;
    }
}
