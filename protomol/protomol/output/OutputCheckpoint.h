/*  -*- c++ -*-  */
#ifndef OUTPUTCHECKPIINT_H
#define OUTPUTCHECKPIINT_H

#include <protomol/output/Output.h>
#include <protomol/base/Timer.h>
#include <string>

namespace ProtoMol {
    class Configuration;

    //____ OutputCheckpoint
    class OutputCheckpoint : public Output {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Constructors, destructors, assignment
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public:
            OutputCheckpoint();
            OutputCheckpoint( const std::string &name, int freq, int start );
        private:
            void WritePositions( int step );
            void WriteVelocities( int step );

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //  From class Output
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        private:
            virtual Output *doMake( const std::vector<Value> &values ) const;
            virtual void doInitialize();
            virtual void doRun( int step );
            virtual void doFinalize( int );
            virtual bool isIdDefined( const Configuration *config ) const;
            virtual bool addDoKeyword() const {return false;}

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // From class Makeable
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public:
            virtual std::string getIdNoAlias() const {return keyword;}
            virtual void getParameters( std::vector<Parameter> & ) const;
            virtual bool adjustWithDefaultParameters( std::vector<Value> &values,
                    const Configuration *config ) const;

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // My data members
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public:
            static const std::string keyword;

        private:
            int mCurrent;
            bool isActive;
            std::string mName;
            std::string mBaseFile;

    };

    template <typename T>
    std::string Append( const std::string& inData, T value ){
        std::ostringstream retStream;

        retStream << inData << value;

        return std::string( retStream.str() );
    }
}

#endif
