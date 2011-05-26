#include <protomol/topology/BuildTopologyFromTpr.h>

#include <protomol/topology/BuildTopology.h>

#include <protomol/base/Exception.h>
#include <protomol/type/String.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/LennardJonesParameterTable.h>
#include <vector>

#include <sstream>
#include <iomanip>

//GROMACS headers etc.
#if defined (HAVE_GROMACS)
extern "C" {

#include "txtdump.h"

}
#endif

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

//use GROMACS exclusions?
#define GROMACSEXCL

// Switch for which Van der Waal radius table to use for GB.
// 0 - Amber default
// 1 - Greg Bowman's modified
#define RADIUS_TABLE 1

//fudge for NETBEANS highlighting
//#define HAVE_GROMACS

//parse functions
bool ProtoMol::parse_iparams(function &func, void * ft, void * ip,
                                ostringstream &os){

    //flag success default
    bool status(true);
    
#if defined (HAVE_GROMACS)

    t_functype *ftype = (t_functype*)ft;

    t_iparams *iparams = (t_iparams*)ip;
    
    switch (*ftype) {
        //bonds
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC:
            func.parameters.push_back(iparams->harmonic.rA);
            func.parameters.push_back(iparams->harmonic.krA);
            func.parameters.push_back(iparams->harmonic.rB);
            func.parameters.push_back(iparams->harmonic.krB);

            os << "b0A= " << setprecision(5) << scientific << setw(12) << iparams->harmonic.rA
                    << ", cbA= " << setw(12) << iparams->harmonic.krA
                    << ", b0B= " << setw(12) << iparams->harmonic.rB
                    << ", cbB= " << setw(12) << iparams->harmonic.krB;

            break;

        //angles
        case F_ANGLES:
        case F_G96ANGLES:
            func.parameters.push_back(iparams->harmonic.rA);
            func.parameters.push_back(iparams->harmonic.krA);
            func.parameters.push_back(iparams->harmonic.rB);
            func.parameters.push_back(iparams->harmonic.krB);

            os << "thA= " << setprecision(5) << scientific << setw(12) << iparams->harmonic.rA
                    << ", ctA= " << setw(12) << iparams->harmonic.krA
                    << ", thB= " << setw(12) << iparams->harmonic.rB
                    << ", ctB= " << setw(12) << iparams->harmonic.krB;

            break;

        //lennard-jones, 1-4
        case F_LJ14:
            func.parameters.push_back(iparams->lj14.c6A);
            func.parameters.push_back(iparams->lj14.c12A);
            func.parameters.push_back(iparams->lj14.c6B);
            func.parameters.push_back(iparams->lj14.c12B);

            os << "c6A= " << setprecision(8) << scientific << setw(15) << iparams->lj14.c6A
                    << ", c12A= " << setw(15) << iparams->lj14.c12A
                    << ", c6B= " << setw(15) << iparams->lj14.c6B
                    << ", c12B= " << setw(15) << iparams->lj14.c12B;

        break;

        //lennard-jones, single
        case F_LJ:
            func.parameters.push_back(iparams->lj.c6);
            func.parameters.push_back(iparams->lj.c12);

            os << "c6A= " << setprecision(8) << scientific << setw(15)
                    << iparams->lj.c6 << ", c12= " << setw(15) << iparams->lj.c12;

            break;

        //RB dihedrals
        case F_RBDIHS:
            //format string
            os << setprecision(8) << scientific;

            //'A' parameters
            for (int i=0; i<NR_RBDIHS; i++){
                func.parameters.push_back(iparams->rbdihs.rbcA[i]);

                if(i!=0) os << ", ";
                os << "rbcA[" << i << "]= " << setw(15) << iparams->rbdihs.rbcA[i];
            }

            //'B' parameters
            for (int i=0; i<NR_RBDIHS; i++){
                func.parameters.push_back(iparams->rbdihs.rbcB[i]);

                os << ", rbcB[" << i << "]= " << setw(15) << iparams->rbdihs.rbcB[i];

            }

            break;

        //'improper' dihedrals
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            func.parameters.push_back(iparams->pdihs.phiA);
            func.parameters.push_back(iparams->pdihs.cpA);
            func.parameters.push_back(iparams->pdihs.phiB);
            func.parameters.push_back(iparams->pdihs.cpB);
            func.parameters.push_back(iparams->pdihs.mult);

            os << setprecision(8) << scientific 
                    << "phiA= " << setw(15) << iparams->pdihs.phiA
                    <<  ", cpA= "<< setw(15) << iparams->pdihs.cpA
                    << ", phiB= " << setw(15) << iparams->pdihs.phiB
                    << ", cpB= " << setw(15) << iparams->pdihs.cpB
                    << ", mult= " << iparams->pdihs.mult;

            break;

        //constraints
        case F_CONSTR:
        case F_CONSTRNC:
            func.parameters.push_back(iparams->constr.dA);
            func.parameters.push_back(iparams->constr.dB);

            os << setprecision(8) << scientific 
                    << "dA= " << setw(15) << iparams->constr.dA
                    << ", dB= " << setw(15) << iparams->constr.dB;
            break;

        //if not found
        default:
            status = false;
            break;
    }

#endif
    
    //return find status
    return status;
}

//main build topology
void ProtoMol::buildTopologyFromTpr(GenericTopology *topo, 
                                    Vector3DBlock &pos, Vector3DBlock &vel,
                                    const string &fname) {
    //define versions of TPR file
    //Version 4.5.3 has tpx_version=73 and includes gb_radius in the tpr file
    enum{ GB_RADII_IN_TPR = 73 };

#if defined (HAVE_GROMACS)

    //----------------------------------------------------------------------
    // Load tpr file
    //----------------------------------------------------------------------

    //load TPR data
    t_state     state;
    rvec        *f = NULL;
    t_inputrec  ir;
    //Tpx struct:
    //int	bIr;		/* Non zero if input_rec is present		*/
    //int	bBox;		/* Non zero if a box is present			*/
    //int	bTop;		/* Non zero if a topology is present		*/
    //int	bX;		/* Non zero if coordinates are present		*/
    //int	bV;		/* Non zero if velocities are present		*/
    //int	bF;		/* Non zero if forces are present		*/
    //int	natoms;		/* The total number of atoms			*/
    //int   ngtc;               /* The number of temperature coupling groups    */
    //real	lambda;		/* Current value of lambda			*/
    t_tpxheader tpx;
    gmx_mtop_t  mtop;

    int file_version, file_generation;

    //load header, in gmxlib/tpxio
    //void read_tpxheader(const char *fn, t_tpxheader *tpx, gmx_bool TopOnlyOK,
    //                int *file_version, int *file_generation)
    read_tpxheader(fname.c_str(),&tpx,TRUE,&file_version,&file_generation);
    report << hint << "File version " << file_version
            << ", generation " << file_generation << "." << endr;

    //load state
    //void read_tpx_state(const char *fn,
    //		    t_inputrec *ir,t_state *state,rvec *f,gmx_mtop_t *mtop)
    read_tpx_state(fname.c_str(),
                 tpx.bIr  ? &ir : NULL,
                 &state,tpx.bF ? f : NULL,
                 tpx.bTop ? &mtop: NULL);

    //test tpx OK
    if ( &tpx == NULL ) {
        THROW("Error loading GROMACS/TPR data.");
    }

    //----------------------------------------------------------------------
    // Load positions and velocities from TPR/Gromacs topology
    //----------------------------------------------------------------------

    //test positions/velocities available
    if( tpx.bX == 0 || tpx.bV == 0 ){
        THROW("No Position or Velocity data.");
    }

    const unsigned int num_atoms = state.natoms;

    //resize the position and velocity arrays
    pos.resize( num_atoms );

    vel.resize( num_atoms );

    //report atom types
    report << debug(810) << "Tpr atom number = " << num_atoms << endr;
    
    //loop over each atom
    for( int i=0; i<num_atoms * 3; i++ ){

        pos.c[i] = state.x[i/3][i%3] * Constant::NM_ANGSTROM; //nm to A; //tpx.bX?
        vel.c[i] = state.v[i/3][i%3] * Constant::NM_ANGSTROM *
                        Constant::TIMEFACTOR * Constant::FS_PS; //nm/ps to A/fs?; //tpx.bV?

    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // First, generate the array of atomtypes
    // Each time a new atom comes up, we need to check if it is
    // already in the vector....
    // NOTE:  this may take a while for large systems; however, it will cut
    // down on the size of the atomTypes vector, and therefore, the amount
    // access time in the back end.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    topo->atoms.clear();
    topo->atomTypes.clear();
    topo->bonds.clear();
    topo->angles.clear();
    topo->dihedrals.clear();
    topo->impropers.clear();

    //Ryckert-Belleman
    topo->rb_dihedrals.clear();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the atom types
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //get atomtypes
    t_atomtypes *atomtypes = &(mtop.atomtypes);
    
    //resize atomtypes
    topo->atomTypes.resize( atomtypes->nr );

    //save size for tests
    const unsigned int atypesize = atomtypes->nr;

    //loop and print data
    for(int i=0;i<atomtypes->nr;i++ ){
        report << debug(810) << "atomtype[" << i << "]={radius=" << atomtypes->radius[i] <<
                ", volume=" << atomtypes->vol[i] <<
                ", gb_radius=" << atomtypes->gb_radius[i] <<
                ", surftens=" << atomtypes->surftens[i] <<
                ", atomnumber=" << atomtypes->atomnumber[i] <<
                ", S_hct=" << atomtypes->S_hct[i] << ")}" << endr;
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the atoms
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // used to get atom names to look up van der waal radii for GB
    vector<string> atomTypeNames;

    // loop over all atoms in the Gromacs topology
    //molecule type topology data
    for(int mt=0; mt<mtop.nmoltype; mt++) {

        gmx_moltype_t *molt = &mtop.moltype[mt];
        gmx_ffparams_t *ffparams = &mtop.ffparams;
        t_atoms *atoms = &(molt->atoms);
        t_atom *atom = atoms->atom;

        //loop over atoms in molecule
        for (int i=0; i<atoms->nr; i++) {

            //get gromacs data
            //atom type, test valid
            const unsigned int atype = atom[i].type;
            if( atype >= atypesize )
                THROW("Undefined atom type.");

            //point to types
            AtomType *tempatomtype = &(topo->atomTypes[atype]);

            //create temporary atom
            Atom tempatom;

            //Update atom type if 'name' not initialized?
            if( tempatomtype->name.length() == 0 ){
                //type not set
                //generate type name from first atom char and index
                stringstream ss;
                ss << atype;
                string str = string(*(atoms->atomname[i])).substr((size_t)0,(size_t)1);
                //
                tempatomtype->name = str + ss.str();
                tempatomtype->mass = atom[i].m;
                tempatomtype->charge = atom[i].q;
                //####just take first char for now.
                tempatomtype->symbolName = str;

                //get radius, test file version
                //if(file_version >= GB_RADII_IN_TPR ){
                    //version 4.5?
                //    tempatomtype->vdwR = atomtypes->gb_radius[atype];
                //}else{
		    /* Removed as part of initial GB fix. 
                    //old version, use orig. params.agb
                    tempatomtype->vdwR = atom_radius( string(*(atoms->atomtype[i])), 1 );

                    if( tempatomtype->vdwR < 0.0 ){
                        THROW("Atom type radius not found.");
                    }
		    */
                //}

            }else{
                //####TODO check new data is consistent? i.e. check new atom mass is the same.
            }

            //report types
            report << debug(810) << "Atom type " << tempatomtype->name << ", " << tempatomtype->mass << ", "
                    << tempatomtype->charge << ", " << tempatomtype->symbolName << ", " << tempatomtype->vdwR << endr;


            // First, we need to find the index. (an integer corresponding
            // to the type of the atom
            tempatom.name = string(*(atoms->atomname[i]));
            tempatom.type = atype;
            tempatom.residue_name = string(*(atoms->resinfo[atom[i].resind].name));
            tempatom.residue_seq = atom[i].resind;
            // Now, the scaled charge.  This is straightforward.
            tempatom.scaledCharge = atom[i].q  * Constant::SQRTCOULOMBCONSTANT;//tempatomtype->charge * Constant::SQRTCOULOMBCONSTANT;
            tempatom.scaledMass = tempatomtype->mass;

	    atomTypeNames.push_back(string(*(atoms->atomtype[i])));
	    
            //report atoms
            report << debug(810) << "Atom " << tempatom.name << ", " <<
                    tempatom.type << ", " << tempatom.residue_name <<
                    ", " << tempatom.residue_seq << ", " << tempatom.scaledCharge <<
                    ", " << tempatom.scaledMass << endr;

            // Now we need the size of the group for heavy atom ordering
            // We need to parse the name for any H's then any numbers following
            // First, if the atom is an H then this is 0
            if (tempatom.name == "H") tempatom.hvyAtom = 0;
            else {
              // Otherwise, we need to parse..
              // Initialize to 1
              tempatom.hvyAtom = 1;
              for (unsigned int pos = 0; pos < tempatom.name.size(); ++pos)
                if (tempatom.name[pos] == 'H') {
                  string number = "";
                  while (isdigit(tempatom.name[++pos]))
                    number += tempatom.name[pos];

                  if (number == "")  // never entered loop, default is 1
                    number = "1";
                  tempatom.hvyAtom += atoi(number.c_str());
                }
            }
            // C/C++ starts at 0, where PSF/PDB at 1
            tempatom.atomNum = i; //####atom->number - 1;
            // Also the molecule - using residue sequence for now
            topo->atoms.push_back(tempatom);

        }

    }

    // calculate the # of degrees of freedom, if there are any bond constraints
    // they will be subtracted later by ModifierShake
    topo->degreesOfFreedom = 3 * topo->atoms.size() - 3;

    report << plain << "D.O.F. = " << topo->degreesOfFreedom << endr;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the functions
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //ffparams
    gmx_ffparams_t * ffparams = &(mtop.ffparams);

    //get fudge parameters
    topo->coulombScalingFactor = ffparams->fudgeQQ;
    topo->LJScalingFactor = 0.5;    //####Actually Read It?

    //vector of functions
    vector<function> functions;

    //loop over all
    for(int i=0; i<ffparams->ntypes; i++){

        //temporary function
        function tempfunc;

        //get name
        tempfunc.name = string( interaction_function[ffparams->functype[i]].name );

        //report
        ostringstream os;
        os << "functype[" << i << "]=" << interaction_function[ffparams->functype[i]].name << ", ";

        //parameters, if successful then save
        if( parse_iparams( tempfunc, (void*)&ffparams->functype[i], 
                            (void*)&ffparams->iparams[i], os ) ){
            //save function
            functions.push_back( tempfunc );
        }else{
            //report << plain << "No matching parms. for " << tempfunc.name << "." << endr;
            THROW("No matching parms. for " + tempfunc.name + ".");
        }

        //end line
        report << debug(800) << os.str() << endr;
        
    }
    
    //report number of functions
    const unsigned int num_functions = functions.size();
    
    report << plain << "Number of functions " << num_functions << "." << endr;

    //force field collections
    vector<ffdata> ffinteract;

    //molecule type topology data
    for(int mt=0; mt<mtop.nmoltype; mt++) {

        for(int j=0; (j<F_NRE); j++) {

            gmx_moltype_t *molt = &mtop.moltype[mt];
            t_functype *functype = ffparams->functype;
            t_ilist *ilist = &molt->ilist[j];
            t_iatom *iatoms;

            //data populated?
            if (ilist->nr > 0) {
                
                report << debug(800) << interaction_function[j].longname << endl << "nr: " << ilist->nr << endr;

                //get pointer to data
                t_iatom *iatoms = ilist->iatoms;

                //loop over data
                for( int i=0; i<ilist->nr; ) {
                    
                    //find type
                    int type = *(iatoms++);
                    //lookup function
                    int ftyp = functype[type];

                    ostringstream os;

                    os << j << " type=" << type << " (" << interaction_function[ftyp].name << ")";

                    //save data
                    ffdata ffd;
                    
                    //name
                    ffd.name = interaction_function[ftyp].name;
                    //long name
                    ffd.longname = interaction_function[j].longname;
                    //function
                    ffd.function = type;//ftyp;

                    //find 'enum' type for testing
                    if( ffd.name.compare("BONDS") == 0 ) ffd.type = BOND;
                    else if( ffd.name.compare("ANGLES") == 0 ) ffd.type = ANGLE;
                    else if( ffd.name.compare("PDIHS") == 0 || ffd.name.compare("PIDIHS") == 0 ) ffd.type = PROPERDIH;
                    else if( ffd.name.compare("RBDIHS") == 0 ) ffd.type = RYCKAERTBELL;
                    else if( ffd.name.compare("LJ14") == 0 ) ffd.type = LJ14;
                    else if( ffd.name.compare("CONSTR") == 0 ) ffd.type = CONSTRAINT;
                    //throw exception and flag if force not identified
                    else{
                        THROW("Unknown GROMAC force type '"+ffd.name+"', '"+ffd.longname+"'.");
                    }

                    //number of atoms associated with force type
                    const unsigned int nratoms = interaction_function[ftyp].nratoms;

                    //atom list
                    for (int k=0; k<nratoms; k++){

                        os << " " << *(iatoms);

                        ffd.atoms.push_back(*(iatoms++));
                    }

                    report << os.str() << endr;

                    //update data counter
                    i += 1 + nratoms;

                    //save interaction
                    ffinteract.push_back(ffd);
                    // "Bond" "Angle" "Proper Dih." "Ryckaert-Bell." "LJ-14"
                }
            }
        }
    }

    //size of force table
    const unsigned int ffsize = ffinteract.size();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the forces
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //This includes the 1-4 exclusions (LJ_14) so creat array here
    //  Resize array.
    vector<exclusion_pair> exclusions;

    // Find the parameters from TPR
    int ignoredBonds = 0;   //preset ignored bonds
    int ignoredAngles = 0;  //and angles

    // get some array sizes
    unsigned int sizeAtomTypes = topo->atomTypes.size();

    //set up LJ table
    topo->lennardJonesParameters.resize(sizeAtomTypes);

    //Setup factors for LJ conversion
    const double tempdistconv6 = pow(Constant::NM_ANGSTROM, 6.0);
    const double C6FACTOR = Constant::KJ_KCAL * tempdistconv6;
    const double C12FACTOR = Constant::KJ_KCAL * tempdistconv6 * tempdistconv6;

    //LJ factor flag
    bool LJfactor_set(false);
    
    //Put LJ data into atomtypes AND fill LJ table with the LJ_SR values
    for( int i=0; i<sizeAtomTypes; i++ ){
        //get function
        const function fnii = functions[i * sizeAtomTypes + i];

        if( fnii.name.compare("LJ_SR") == 0 ){

            const double c6 = fnii.parameters[0];
            const double c12 = fnii.parameters[1];

            if( c6 != 0.0 && c12 != 0.0 ){
                //find sigma frpm c6/c12
                const double sigma = pow(c12/c6, (1.0/6.0)) * Constant::NM_ANGSTROM;
                topo->atomTypes[i].sigma = sigma;
 
                //find epsilon from c6/c12
                topo->atomTypes[i].epsilon = c6*c6/(4.0*c12) * Constant::KJ_KCAL;
                
            }else{
                report << plain << "Warning: LJ_SR Function has zero c6 or c12 parameter, setting sigma=1 and epsilon=0." << endr;
                topo->atomTypes[i].sigma = 1.0;
                topo->atomTypes[i].epsilon = 0.0;
            }

            //copy to 1-4
            topo->atomTypes[i].sigma14 = topo->atomTypes[i].sigma;
            topo->atomTypes[i].epsilon14 = topo->atomTypes[i].epsilon;
            
        }else{
            THROW("Error: LJ_SR Function is of type '"+fnii.name+"'.");
        }

        //Now fill table
        for( int j=0; j<sizeAtomTypes; j++ ){

            //get function
            const function fn = functions[i * sizeAtomTypes + j];

            if( fn.name.compare("LJ_SR") == 0 ){

                //parameters
                LennardJonesParameters paramsij;

                paramsij.A = fn.parameters[1] * C12FACTOR; //c12
                paramsij.B = fn.parameters[0] * C6FACTOR; //c6
                paramsij.A14 = paramsij.A * topo->LJScalingFactor;
                paramsij.B14 = paramsij.B * topo->LJScalingFactor;

                topo->lennardJonesParameters.set(i, j, paramsij);

            }else{
                THROW("Error: LJ_SR Function is of type '"+fn.name+"'.");
            }

        }
    }

    //loop over interaction list
    for( int ffii=0; ffii < ffsize; ffii++ ){
        
        //get data
        const ffdata ffd = ffinteract[ffii];
        
        //get function
        const function fn = functions[ffd.function];

        //switch on each force type
        switch(ffd.type){

            //~~~~Bonds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case BOND: {
                    // if we have found this bond type then copy the bond parameters
                    // into the topology

                    //test data sizes match
                    if( fn.parameters.size() < 2 || ffd.atoms.size() != 2 )
                            THROW("Error in BOND function '"+fn.name+"', '"+ffd.name+"'.");

                    Bond tempbond;
                    tempbond.restLength = fn.parameters[0] * Constant::NM_ANGSTROM;
                    tempbond.springConstant = fn.parameters[1]
                                                * Constant::KJ_KCAL * Constant::ANGSTROM_NM * Constant::ANGSTROM_NM * 0.5;
                    tempbond.atom1 = ffd.atoms[0];
                    tempbond.atom2 = ffd.atoms[1];
                    topo->bonds.push_back(tempbond);

                    // populate the vector of bonds maintained at each atom
                    topo->atoms[tempbond.atom1].mybonds.push_back((topo->bonds.size()) - 1);
                    topo->atoms[tempbond.atom2].mybonds.push_back((topo->bonds.size()) - 1);

                    if (tempbond.springConstant == 0.0) ignoredBonds++;
                }
                //
                break;

            //~~~~Angles~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case ANGLE: {
                    // if we have found this angle type then copy the angle parameters
                    // into the topology

                    //test data sizes match
                    if( fn.parameters.size() < 2 || ffd.atoms.size() != 3 )
                            THROW("Error in ANGLE function '"+fn.name+"', '"+ffd.name+"'.");

                    Angle tempangle;
                    //####note fliping of atoms
                    tempangle.atom1 = ffd.atoms[0];
                    tempangle.atom2 = ffd.atoms[1];
                    tempangle.atom3 = ffd.atoms[2];
                    tempangle.restAngle = fn.parameters[0] * M_PI/180.0;
                    tempangle.forceConstant = fn.parameters[1]
                                                * Constant::KJ_KCAL * 0.5; //times 1/2 as Amber is 1/2 k(a-a_0)^2;
                    // no Urey-Bradley term specified
                    tempangle.ureyBradleyConstant = 0.0;
                    tempangle.ureyBradleyRestLength = 0.0;
                    topo->angles.push_back(tempangle);
                    if(tempangle.forceConstant == 0.0) ignoredAngles++;
                }
                //
                break;

            //~~~~Proper Dihedrals~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case PROPERDIH: {
                    // if we have found this proper dihedral type then copy the
                    // proper dihedral parameters into the topology
                
                    //test data sizes match
                    if( fn.parameters.size() < 5 || ffd.atoms.size() != 4 )
                            THROW("Error in PROPERDIH function '"+fn.name+"', '"+ffd.name+"'.");

                    Torsion torsion;
                    torsion.atom1 = ffd.atoms[0];
                    torsion.atom2 = ffd.atoms[1];
                    torsion.atom3 = ffd.atoms[2];
                    torsion.atom4 = ffd.atoms[3];
                    torsion.periodicity.push_back((int)fn.parameters[4]); //equiv to mult.
                    torsion.phaseShift.push_back(fn.parameters[0] * M_PI / 180.0); //phiA
                    torsion.forceConstant.push_back(fn.parameters[1] * Constant::KJ_KCAL); //cpA
                    torsion.multiplicity = 1;
                    topo->dihedrals.push_back(torsion);                
                }
                //
                break;

            //~~~~Ryckaert-Bell Dihedrals~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case RYCKAERTBELL: {
                    // if we have found this RB diheral type then copy the
                    // dihedral parameters into the topology

                    //test data sizes match
                    if( fn.parameters.size() < 6 || ffd.atoms.size() != 4 )
                            THROW("Error in RYCKAERTBELL function '"+fn.name+"', '"+ffd.name+"'.");

                    // if we have found this dihedral type then copy the
                    // dihedral parameters into the topology
                    RBTorsion rb_torsion;
                    rb_torsion.atom1 = ffd.atoms[0];
                    rb_torsion.atom2 = ffd.atoms[1];
                    rb_torsion.atom3 = ffd.atoms[2];
                    rb_torsion.atom4 = ffd.atoms[3];

                    rb_torsion.C0  = fn.parameters[0] * Constant::KJ_KCAL;
                    rb_torsion.C1  = fn.parameters[1] * Constant::KJ_KCAL;
                    rb_torsion.C2  = fn.parameters[2] * Constant::KJ_KCAL;
                    rb_torsion.C3  = fn.parameters[3] * Constant::KJ_KCAL;
                    rb_torsion.C4  = fn.parameters[4] * Constant::KJ_KCAL;
                    rb_torsion.C5  = fn.parameters[5] * Constant::KJ_KCAL;

                    topo->rb_dihedrals.push_back(rb_torsion);
                }
                //
                break;

            //~~~~Lennard-Jones~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //over lay with LJSR? in A? and A14?
            case LJ14: {
                    // if we have found this LJ type then copy the
                    // LJ parameters into the topology

                    //save pair
                    exclusions.push_back(exclusion_pair( ffd.atoms[0], ffd.atoms[1]) );
                    
                    //test data sizes match
                    if( fn.parameters.size() != 4 || ffd.atoms.size() != 2 )
                            THROW("Error in LJ14 function '"+fn.name+"', '"+ffd.name+"'.");

                    //parameters
                    LennardJonesParameters paramsij;
                    
                    paramsij.A14 = fn.parameters[1] * C12FACTOR; //c12
                    paramsij.B14 = fn.parameters[0] * C6FACTOR; //c6

                    //save in table, get types first
                    const int atype1 = topo->atoms[ffd.atoms[0]].type;
                    const int atype2 = topo->atoms[ffd.atoms[1]].type;
                    
                    //get non-exclusion functions
                    const function fnNonExcl = functions[atype1 * sizeAtomTypes + atype2];

                    paramsij.A = fnNonExcl.parameters[1] * C12FACTOR; //c12
                    paramsij.B = fnNonExcl.parameters[0] * C6FACTOR; //c6

                    if( !LJfactor_set ){
                        LJfactor_set = true;
                        topo->LJScalingFactor = paramsij.A14 / paramsij.A;

                        report << debug(800) << "LJ Ratio set to " << topo->LJScalingFactor << "." << endr;

                    }
                    
                    //save
                    topo->lennardJonesParameters.set(atype1, atype2,
                                                        paramsij);
                }
                //
                break;

        }
    }

    //ignored bonds etc.?
    if (ignoredBonds > 0)
        report << hint << "Systems contains " << ignoredBonds
               << " bonds with zero force constants." << endr;

    if (ignoredAngles > 0)
        report << hint << "Systems contains " << ignoredAngles <<
            " angles with zero force constants." << endr;

    //check LJ table
    for( int i=0; i<sizeAtomTypes; i++ ){
        for( int j=i+1; j<sizeAtomTypes; j++ ){

            double sigma_i = topo->atomTypes[i].sigma;
            double sigma_j = topo->atomTypes[j].sigma;
            double sigma14_i = topo->atomTypes[i].sigma14;
            double sigma14_j = topo->atomTypes[j].sigma14;

            double epsilon_i = topo->atomTypes[i].epsilon;
            double epsilon_j = topo->atomTypes[j].epsilon;
            double epsilon14_i = topo->atomTypes[i].epsilon14;
            double epsilon14_j = topo->atomTypes[j].epsilon14;

            double r_ij, e_ij, r14_ij, e14_ij;
            r_ij = 0.5 * (sigma_i + sigma_j);
            e_ij = sqrt(epsilon_i * epsilon_j);
            r14_ij = 0.5 * sigma14_i + sigma14_j;
            e14_ij = sqrt(epsilon14_i * epsilon14_j);

            const LennardJonesParameters paramsij = topo->lennardJonesParameters(i,j);

            ostringstream os;

            os << "LJ table A " << i << " " << j << " " << paramsij.A << " " << power<12>(r_ij) * e_ij * 4.0;
            os << ", LJ table B " << i << " " << j << " " << paramsij.B << " " << power<6>(r_ij) * e_ij * 4.0;
            ///#### FudgeLJ=0.5 should be read from the parameter files
            os << ", LJ table A14 " << i << " " << j << " " <<  paramsij.A14 << " " << topo->LJScalingFactor * paramsij.A;
            os << ", LJ table B14 " << i << " " << j << " " <<  paramsij.B14 << " " << topo->LJScalingFactor * paramsij.B << endl;

            report << debug(800) << os.str();
        }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the exclusions
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //creat vector of exclusions
    vector<exclusion_pair> all_exclusions;

    //  calc. exclusion size
    unsigned int exclsize = 0;

    //find each molecule exclusion size
    for( int mt=0; mt<mtop.nmoltype; mt++){

        gmx_moltype_t *molt = &mtop.moltype[mt];

        //point to blocks
        t_blocka *block = &molt->excls;

        exclsize += ( block->nra - block->nr ) / 2;
    }

    //do resize
    //##topo->exclusions.resize(exclsize);

    //Now get the exclusions, over each molecule
    for( int mt=0; mt<mtop.nmoltype; mt++){

        gmx_moltype_t *molt = &mtop.moltype[mt];

        //point to blocks
        t_blocka *block = &molt->excls;

        //control variables
        int start, end;

        end = start = 0;

        //check start
        if( block->index[start] != 0 ){
            THROW("block->index[0] should be 0");
        }else{
            //find exclusions
            for (int i=0; i<block->nr; i++){
                end=block->index[i+1];
                
                //report?
                ostringstream os;

                if (end<=start){

                    os << "excls[" << i << "]={";

                }else{

                    os << "excls[" << i << "][" << start << ".." << end << "]={";

                }
                for (int j=start; j<end; j++){

                    if (j>start) os << ", ";
                    os << block->a[j];
                    
                    //add exclusion?
                    const int secondatom = block->a[j];

                    //save unique only
                    if( secondatom > i ){
                        //save pair
                        all_exclusions.push_back( exclusion_pair( i, secondatom ) );

                        //##topo->exclusions.add(i, block->a[j], EXCLUSION_FULL);
                    }

                }

                report << debug(800) << os.str() << "}" << endr;

                start = end;
            }
        }

        //report exclusions
        report << debug(800) << "Exclusion size: " << all_exclusions.size() << "(" << exclsize << ")."<< endr;
        
        //check correct termination
        if( end != block->nra ){
            THROW("Exclusion tables inconsistent.");
        }
        
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Generalized Born data
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //test for flag
    if ( ir.implicit_solvent == 1 ) {
        //set parameters from topology
        topo->doGBSAOpenMM = 1;
        topo->implicitSolvent = GBSA;

        topo->obcType = ir.gb_algorithm + 1;
        topo->alphaObc = ir.gb_obc_alpha;
        topo->betaObc = ir.gb_obc_beta;
        topo->gammaObc = ir.gb_obc_gamma;
        topo->dielecOffset = ir.gb_dielectric_offset * Constant::NM_ANGSTROM;

 
        report << debug(800) << "Implicit solvent " << ir.implicit_solvent
                << ", OBC type " << topo->obcType
                << ", alpha " << topo->alphaObc
                << ", beta " << topo->betaObc
                << ", gamma " << topo->gammaObc
                << ", dielec offset " << topo->dielecOffset
                << "." << endr;
       
        //setup atom structs
        unsigned int atomsSize = topo->atoms.size();
        //string hname("H"), cname("C"),
        for (unsigned int i=0;i < atomsSize; i++) {

            Atom *tempatom = &(topo->atoms[i]);

            tempatom->myGBSA_T = new GBSAAtomParameters();

            int type = tempatom->type;

            string name = topo->atomTypes[type].name;

	    tempatom->myGBSA_T->vanDerWaalRadius = atom_radius( atomTypeNames[i], RADIUS_TABLE);

	    if(tempatom->myGBSA_T->vanDerWaalRadius < 0.0) {
		    stringstream ss;
		    ss << "Radius not found for Atom type " << atomTypeNames[i];
		    THROW(ss.str().c_str());
	    }


            tempatom->myGBSA_T->offsetRadius = 0.09;
            //tempatom->myGBSA_T->scalingFactor
            if ( name[0] == 'H') {
               tempatom->myGBSA_T->scalingFactor = 0.85;
            }else if ( name[0] == 'C') {
               tempatom->myGBSA_T->scalingFactor = 0.72;
            }else if ( name[0] == 'N') {
               tempatom->myGBSA_T->scalingFactor = 0.79;

            }else if ( name[0] == 'O') {
               tempatom->myGBSA_T->scalingFactor = 0.85;
            } else if ( name[0] == 'P') {
               tempatom->myGBSA_T->scalingFactor = 0.86;
            } else if (name[0] == 'S') {
               tempatom->myGBSA_T->scalingFactor = 0.96;
            } else {
               tempatom->myGBSA_T->scalingFactor = 0.8;
            }

            //allocate the array to store derivatives of born radius w.r.t. r_{ij}'s
            if (tempatom->myGBSA_T->bornRadiusDerivatives == NULL) {
                tempatom->myGBSA_T->SetSpaceForBornRadiusDerivatives(atomsSize);
            }

            if (tempatom->myGBSA_T->Lvalues == NULL) {
                tempatom->myGBSA_T->SetSpaceLvalues(atomsSize);
            }
            if (tempatom->myGBSA_T->Uvalues == NULL) {
                tempatom->myGBSA_T->SetSpaceUvalues(atomsSize);
            }

            if (tempatom->myGBSA_T->distij == NULL) {
               tempatom->myGBSA_T->SetSpaceDistij(atomsSize);
            }
        }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Tables and exclusions
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // store the molecule information
    buildMoleculeTable(topo);

    //Calculate exclusions based on TPR file rather than physical topology?
#if defined(GROMACSEXCL)
    
    //  Resize array to atoms size.
    topo->exclusions.resize(topo->atoms.size());

    //add any 1-4 defined in the TPR file
    const unsigned int excl_size = exclusions.size();

    for( int i=0; i<excl_size; i++ ){
        //set
        topo->exclusions.add(exclusions[i].atom1, exclusions[i].atom2, EXCLUSION_MODIFIED);
    }

    //and full exclusions
    const unsigned int all_exclsize = all_exclusions.size();
        
    for( int i=0; i<all_exclsize; i++ ){
        if( topo->exclusions.check( all_exclusions[i].atom1, all_exclusions[i].atom2 ) != EXCLUSION_MODIFIED )
            topo->exclusions.add(all_exclusions[i].atom1, all_exclusions[i].atom2,
                EXCLUSION_FULL);
    }
    
#else
    //build exclusions (not using exclusions from above, recalc from topology)
    buildExclusionTable(topo, topo->exclude);


    //add any defined in the TPR file not associated with a dihedral
    const unsigned int excl_size = exclusions.size();

    for( int i=0; i<excl_size; i++ ){
        //check if already set
        if( topo->exclusions.check( exclusions[i].atom1, exclusions[i].atom2 ) != EXCLUSION_MODIFIED )
            topo->exclusions.add(exclusions[i].atom1, exclusions[i].atom2, EXCLUSION_MODIFIED);
    }
#endif

    //optimize again
    topo->exclusions.optimize();
  
#endif
  
}

//atomic radius from lookup
//##TODO:Should be read from 'params.agb'
double ProtoMol::atom_radius( string atom_type, int set ){

    enum{ ARRAYSZ = 96};
    
string amber_atom[2][ARRAYSZ] =
{
 //Original set
 {"amber99_0", "amber99_2", "amber99_3", "amber99_4", "amber99_5",
 "amber99_6", "amber99_7", "amber99_8", "amber99_9", "amber99_10",
 "amber99_11", "amber99_12", "amber99_13", "amber99_14", "amber99_17",
 "amber99_18", "amber99_19", "amber99_20", "amber99_21", "amber99_22",
 "amber99_23", "amber99_24", "amber99_25", "amber99_26", "amber99_27",
 "amber99_28", "amber99_34", "amber99_35", "amber99_36", "amber99_37",
 "amber99_38", "amber99_39", "amber99_40", "amber99_41", "amber99_42",
 "amber99_43", "amber99_44", "amber99_45", "amber99_47", "amber94_0",
 "amber94_1", "amber94_2", "amber94_3", "amber94_4", "amber94_5",
 "amber94_6", "amber94_7", "amber94_8", "amber94_9", "amber94_10",
 "amber94_11", "amber94_12", "amber94_13", "amber94_14", "amber94_15",
 "amber94_16", "amber94_17", "amber94_17b", "amber94_18", "amber94_19",
 "amber94_20", "amber94_21", "amber94_22", "amber94_23", "amber94_24",
 "amber94_25", "amber94_26", "amber94_27", "amber94_28", "amber94_29",
 "amber94_30", "amber94_31", "amber94_33", "amber94_34", "amber94_35",
 "amber94_36", "amber94_37", "amber94_38", "amber94_39", "amber94_40",
 "amber94_41", "amber94_42", "amber94_43", "amber94_44", "amber94_45",
 "amber94_46", "amber94_47", "amber94_48", "amber94_54", "amber94_55",
 "amber94_60", "amber94_61", "amber94_62", "amber94_63", "amber94_64",
 "amber94_65" },
 //Greg Bowman set
 {"amber99_0","amber99_2","amber99_3","amber99_4","amber99_5","amber99_6",
 "amber99_7","amber99_8","amber99_9","amber99_10","amber99_11","amber99_12",
 "amber99_13","amber99_14","amber99_17","amber99_18","amber99_19",
 "amber99_20","amber99_21","amber99_22","amber99_23","amber99_24",
 "amber99_25","amber99_26","amber99_27","amber99_28","amber99_34",
 "amber99_35","amber99_36","amber99_37","amber99_38","amber99_39",
 "amber99_40","amber99_41","amber99_42","amber99_43","amber99_44",
 "amber99_45","amber99_47","amber94_0","amber94_2","amber94_3","amber94_4",
 "amber94_5","amber94_6","amber94_7","amber94_8","amber94_9","amber94_10",
 "amber94_11","amber94_12","amber94_13","amber94_14","amber94_17",
 "amber94_18","amber94_19","amber94_20","amber94_21","amber94_22",
 "amber94_23","amber94_24","amber94_25","amber94_26","amber94_27",
 "amber94_28","amber94_34","amber94_35","amber94_36","amber94_37",
 "amber94_38","amber94_39","amber94_40","amber94_41","amber94_42",
 "amber94_43","amber94_44","amber94_45","amber94_47",
 "","","","","","","","","","","","","","","","","","" }
};

double rad[2][ARRAYSZ] =
{
 //Original set
 {1.875, 1.875, 1.875, 1.875, 1.875, 1.875, 1.875, 1.875,
 1.875, 1.875, 1.900, 1.875, 1.875, 1.875, 1.150, 1.250,
 1.250, 1.250, 1.250, 1.250, 1.250, 1.250, 1.050, 1.250,
 1.050, 1.250, 1.706, 1.706, 1.706, 1.706, 1.706, 1.625,
 1.706, 1.480, 1.535, 1.535, 1.535, 1.480, 1.775, 1.200,
 1.500, 1.700, 1.700, 1.700, 1.700, 1.700, 1.700, 1.700,
 1.700, 1.700, 1.700, 1.700, 1.700, 1.700, 1.500, 1.500,
 1.300, 1.300, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200,
 1.200, 1.200, 1.200, 1.200, 1.200, 1.500, 1.700, 1.500,
 1.500, 1.550, 1.550, 1.550, 1.550, 1.550, 1.550, 1.550,
 1.500, 1.500, 1.500, 1.500, 1.500, 1.850, 1.800, 1.800,
 1.500, 1.200, 1.200, 1.500, 1.200, 1.500, 1.200, 1.500},
 //Greg Bowman set
 {1.875,1.875,1.875,1.875,1.875,1.875,1.875,1.875,1.875,
 1.875,1.900,1.875,1.875,1.875,1.150,1.250,1.250,1.250,
 1.250,1.250,1.250,1.250,1.050,1.250,1.050,1.250,1.706,
 1.706,1.706,1.706,1.706,1.625,1.706,1.480,1.535,1.535,
 1.535,1.480,1.775,1.875,1.875,1.875,1.875,1.875,1.875,
 1.875,1.875,1.875,1.875,1.900,1.875,1.875,1.875,1.150,
 1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.050,1.250,
 1.050,1.250,1.706,1.706,1.706,1.706,1.706,1.625,1.706,
 1.480,1.535,1.535,1.535,1.480,1.775,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
};

    //search array
    for( int i=0; i<ARRAYSZ; i++ ){
        if( atom_type.compare(amber_atom[set][i]) == 0 ){
            return rad[set][i];
        }
    }

    //flag error if not in list
    return -1.0;

}
