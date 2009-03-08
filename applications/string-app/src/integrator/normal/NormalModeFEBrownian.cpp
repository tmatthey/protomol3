#include <src/integrator/normal/NormalModeFEBrownian.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#include <src/integrator/normal/StringModifierForceProjection.h>


using namespace std;

using namespace ProtoMol::Report;

using std::string;
using std::vector;


namespace ProtoMol {
	//__________________________________________________ NormalModeFEBrownian

	const string NormalModeFEBrownian::keyword( "NormalModeFEBrownian" );

	NormalModeFEBrownian::NormalModeFEBrownian() : STSIntegrator(), StringNormalModeUtilities() 
	{
		myWriter = NULL;
		myWriter2 = NULL;	

	}

	NormalModeFEBrownian::NormalModeFEBrownian(Real timestep, int firstmode, int nummode, Real gamma, int seed, Real temperature,  
		std::string avff, std::string inff,//####added avff, inff for diagnostics
		ForceGroup *overloadedForces) 
		: STSIntegrator(timestep,overloadedForces), StringNormalModeUtilities( firstmode, nummode, gamma, seed, temperature), 
		avForceFile(avff), inForceFile(inff) //####added avForceFile, inForceFile for diagnostics
	{
		myWriter = NULL;
		myWriter2 = NULL;
	}

	NormalModeFEBrownian::~NormalModeFEBrownian() 
	{  
		if(myWriter != NULL) delete myWriter;
		if(myWriter2 != NULL) delete myWriter2;
	}

	void NormalModeFEBrownian::initialize(ProtoMolApp* app){
			STSIntegrator::initialize(app);
			initializeForces();
			//NM initialization if OK
			StringNormalModeUtilities::initialize((int)app->positions.size(), app->topology, myForces, NO_NM_FLAGS); //last for non-complimentary forces
			//
			//initialize minimizer noise vars
			randStp = 0.0;
			//zero instantaneous and average force Vector3DBlock
			tempV3DBlk.resize(_N);
			temp2V3DBlk.resize(_N);
			//***********************************************************	
			//####diagnostics
			if(avForceFile != ""){
				myWriter = new XYZTrajectoryWriter();
				if(!myWriter->openWith(avForceFile, app->topology->atoms, app->topology->atomTypes))
					report << error << "Can't open output file '"<<avForceFile<<"'."<<endr;
			}

			if(inForceFile != ""){
				myWriter2 = new XYZTrajectoryWriter();
				if(!myWriter2->openWith(inForceFile, app->topology->atoms, app->topology->atomTypes))
					report << error << "Can't open output file '"<<inForceFile<<"'."<<endr;
			}

	}

	void NormalModeFEBrownian::run(int numTimesteps) {
		Real h = getTimestep() * Constant::INV_TIMEFACTOR;
		Real actTime;

		if( numTimesteps < 1 )
			return;

		//check valid eigenvectors
		if(*Q == NULL)
			report << error << "No Eigenvectors for NormalMode integrator."<<endr;
		//time calculated in forces! so fix here
		actTime = app->topology->time + numTimesteps * getTimestep();
		//main loop
		//zero average force
		tempV3DBlk.zero(_N);
		aveForceCount = 0;
		//
		for( int i = 0; i < numTimesteps; i++ ) {
                        //cout<<"******START MAIN LOOP OF NMBROWNIAN"<<endl;
			//****main loop*************************************
			//
			preStepModify();
			//      
			calculateForces();
                        //for( unsigned int j = 0; j < app->positions.size(); ++j) cout<<"In Brownian : atom "<<j<<" "<<(*myForces)[j]<<endl;
                        //cout<<"PE from NormalModeFEBrownian::run "<<app->energies.potentialEnergy()<<endl;
			//
			for( unsigned int j = 0; j < app->positions.size(); ++j)
				app->positions[j]  += (*myForces)[j] * h / (app->topology->atoms[j].scaledMass * myGamma);
			//add random force
			genProjGauss(&gaussRandCoord1, app->topology);
			randStp = sqrt(2 * Constant::BOLTZMANN * myTemp * h / myGamma);
                        //cout<<"NormalModeFEBrownian : randStp "<<randStp<<endl;
			app->positions.intoWeightedAdd(randStp,gaussRandCoord1);
			//
			//
			postStepModify();
                        //cout<<"******STOP MAIN LOOP OF NMBROWNIAN"<<endl;
		}	
		//fix average, and output
		if(aveForceCount){
                        //cout<<"Brownian : aveForceCount "<<aveForceCount<<endl;
			for( unsigned int i=0;i < app->positions.size(); i++) tempV3DBlk[i] /= (Real)aveForceCount;
                        //for( unsigned int i=0;i < app->positions.size(); i++) cout<<"IN Brownian : mean force for atom "<<i<<" "<<tempV3DBlk[i]<<endl;
			//####diagnostics
			if(avForceFile != "")
				myWriter->write(tempV3DBlk);

			//####

		}
		//fix time
		app->topology->time = actTime;
		//
	}  

	void NormalModeFEBrownian::getParameters(vector<Parameter>& parameters) const {
		STSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NoConstraints()),-1,Text("First mode to use in set")));
		parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NoConstraints()),-1,Text("Number of modes propagated")));
		parameters.push_back(Parameter("gamma",Value(myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
		parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234,Text("Langevin random seed")));
		parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
		//####diagnostics
		parameters.push_back(Parameter("avForceFile",Value(avForceFile,ConstraintValueType::NoConstraints()),std::string(""),Text("Average force filename")));
		parameters.push_back(Parameter("inForceFile",Value(inForceFile,ConstraintValueType::NoConstraints()),std::string(""),Text("Instantaneous force filename")));

		//####
	}

	STSIntegrator* NormalModeFEBrownian::doMake(const vector<Value>& values,ForceGroup* fg)const{
		return new NormalModeFEBrownian(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7], //####last 2 for diagnostics
		fg);
	}

	void NormalModeFEBrownian::addModifierAfterInitialize(){
		adoptPostForceModifier(new StringModifierForceProjection(this));
		STSIntegrator::addModifierAfterInitialize();
	}

	//override force projection for post force modifier, and find average in complement space
	void NormalModeFEBrownian::forceProjection(){
		unsigned int count = myForces->size();
		if((*Q) != NULL){
			// myForces has total forces
			for( unsigned int i=0;i < count; i++) temp2V3DBlk[i] = (*myForces)[i];
			// project myForces onto fast subspace 
			subspaceForce(myForces, myForces);
			// difference between old myForces stored in temp2V3DBlk and myForces
			// gives us instantaneous slow force
			for( unsigned int i=0;i < count; i++) temp2V3DBlk[i] -= (*myForces)[i];
                        //for( unsigned int j = 0; j < app->positions.size(); ++j) cout<<"In Brownian in force projection [22]: atom "<<j<<" "<<temp2V3DBlk[j]<<endl;
			//###diagnostics
			if(inForceFile != "")
				myWriter2->write(temp2V3DBlk);

			// add slow force to running sum in tempV3DBlk[i] 
			for( unsigned int i=0;i < count; i++) tempV3DBlk[i] += temp2V3DBlk[i];
                        //tempV3DBlk.intoAdd(temp2V3DBlk);
			aveForceCount++;
                        //for( unsigned int j = 0; j < app->positions.size(); ++j) cout<<"In Brownian in force projection (accumulated mean force) : atom "<<j<<" "<<tempV3DBlk[j]<<endl;

		}
               
	}


}

