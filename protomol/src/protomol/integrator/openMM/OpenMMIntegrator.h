#ifndef OPENMMINTEGRATOR_H
#define OPENMMINTEGRATOR_H

#include <protomol/integrator/STSIntegrator.h>

#include <OpenMM.h>

#ifdef HAVE_OPENMM_LTMD
#include <LTMD/Integrator.h>
#include <LTMD/Parameters.h>
#endif

#ifdef HAVE_OPENMM_FBM
#include <OpenMMFBM/FBMParameters.h>
#endif

namespace ProtoMol {
	class ScalarStructure;
	class ForceGroup;
	class Vector3DBlock;

	class OpenMMIntegrator : public STSIntegrator {
		public:
			OpenMMIntegrator();
			OpenMMIntegrator( const std::vector<Value>& params, ForceGroup *overloadedForces );
			~OpenMMIntegrator();
		public:
			virtual std::string getIdNoAlias() const { return keyword; }
			virtual void getParameters( std::vector<Parameter> &parameters ) const;
			virtual unsigned int getParameterSize() const;

			virtual void initialize( ProtoMolApp *app );
			virtual long run( const long numTimesteps );
		private:
			virtual STSIntegrator *doMake( const std::vector<Value> &values, ForceGroup *fg ) const;
			void getObcScaleFactors( std::vector<double>& scaleFactors );
		public:
			static const std::string keyword;
		protected:
			// OpenMM Parameters
			int mPlatform, mPropagatorDevice, mBlockDevice, mMinSteps;
			Real mTolerance;

			OpenMM::System *system;
			OpenMM::Context *context;
			OpenMM::Integrator *integrator;

			// Integrator Parameters
			int mSeed;
			Real mTemperature, mGamma;

			// Solvent Parameters
			int mCommonMotionRate;
			Real mGBSAEpsilon, mGBSASolvent;

			// CutoffParams
		        Real mNonbondedCutoff, mGBCutoff;

			// Force Switches
			bool isUsingHarmonicBondForce;
			bool isUsingHarmonicAngleForce;
			bool isUsingNonBondedForce;
			bool isUsingRBDihedralForce;
			bool isUsingPeriodicTorsionForce;
			bool isUsingGBForce;
			bool isUsingSCPISMForce;
			bool isUsingImproperTorsionForce;
			bool isUsingUreyBradleyForce;
			bool isUsingCHARMMLennardJonesForce;

			// LTMD Data
			bool isLTMD;
			std::vector<std::string> mForceList;

#ifdef HAVE_OPENMM_LTMD
			OpenMM::LTMD::Parameters mLTMDParameters;
#endif

#ifdef HAVE_OPENMM_FBM
			FBMParameters mFBMParameters;
#endif

	};
}

#endif // OPENMMINTEGRATOR_H
