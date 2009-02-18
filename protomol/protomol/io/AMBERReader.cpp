#include <protomol/io/AMBERReader.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>

using std::string;
using std::vector;
using std::endl;
using std::cout;
using std::find;
using std::stringstream;
using namespace ProtoMol::Report;

//#define DEBUG_AMBER

namespace ProtoMol {

  //_________________________________________________________________AMBERReader

  AMBERReader::AMBERReader():     
    Reader(),myPSF(NULL),myPAR(NULL){}

  AMBERReader::AMBERReader(const std::string& filename):
    Reader(filename),myPSF(NULL),myPAR(NULL){} 

  AMBERReader::~AMBERReader(){
    if(myPSF != NULL)
      delete myPSF;
    if(myPAR != NULL)
      delete myPAR;
    
  }

  //amber topology does not start with any keyword as in psfs.
  //That's why there is nothing to do to check format.	
  bool AMBERReader::tryFormat(){    //tryFormat check file format
    if(!open())
      return false;
    else return true;
    //string topHead;
    //file >> topHead;
    //return (file.good() && equalNocase("AMBER",topHead));
  }

  bool AMBERReader::read() {
    std::cout<<filename<<std::endl;
    if(myPSF == NULL) {
      myPSF = new PSF();
      std::cout<<"Creating PSF structure"<<endl;}	
    if(myPAR == NULL){
      myPAR = new PAR();
      std::cout<<"Creating PAR structure"<<endl;}
    return read(*myPSF, *myPAR);
  }

  bool AMBERReader::read(PSF& psf,PAR& par) {
    
    //check if the file is open
    if ( !open()) {
      cout<<"Not open file:false"<<endl;
      return false;
    }

    //The file is open... so we can start reading topology information.
    //first clear the data structures.
    psf.clear();
    par.clear();
    

    //Now proceed and read the toplogy information	
    string line = getline();
    //cout<<line<<endl;
    vector<int> ptrVals;
    int numRecords, numAngles;
    
    stringstream ss;	
    string str;
    string keyword;
    vector<int> iac;
    vector<int> ico;
    vector<string> residues;
    vector<int> residue_pointers;
    vector<Real> cn1;
    vector<Real> cn2;
    vector<Real> rk;
    vector<Real> req;
    vector<Real> tk;
    vector<Real> teq;
    vector<Real> phi_k;
    vector<Real> phi_n;
    vector<Real> phi_phase;
    vector<Real> hbond_acoff;
    vector<Real> hbond_bcoff;

    line = getline();
    //cout<<line<<endl;

    while(!file.eof()) {
      //ignore blank lines first.
      int temp_atom_number;

      //line =getline();
      //cout<<"Inside while loop "<<line<<endl;

      if (isBlank(line)) {
    line = getline();
    continue; 
      }
      else {
    //extract the header and check
    if(std::find(line.begin(),line.end(),'%') !=line.begin()) return false;
    else {	
      stringstream ss(line);
      ss >> keyword;
      if (equal("%FLAG", keyword)) {ss >> keyword; 
      //std::cout<<"Keyword = "<<keyword<<std::endl;
      }
      else if(std::find(line.begin(),line.end(),'%') !=line.begin()) return false; 
      line = getline();
    }	
      }
      if(equalStartNocase("POINTERS",keyword)) {
    //line = getline();
    while((!line.empty()) && (std::find(line.begin(), line.end(),'%') == line.begin())) {
      line = getline();
      //cout<<line<<endl;
    }
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      line = removeBeginEndBlanks(line);
      ss.clear();
      ss << line;
      //std::cout<<line<<std::endl;
      while(ss>>str) {
        //cout<<str<<"\t";	
        ptrVals.push_back(toInt(str));
      }
      //cout<<ptrVals.size()<<endl;
      line = getline();
    }
    //cout<<ptrVals.size()<<endl;
    //cout<<ptrVals[4]<<endl;
    //cout<<line<<endl;
      }
      else if(equalStartNocase("ATOM_NAME",keyword)){

    //line = getline();
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    numRecords = ptrVals[0];
    //cout<<"NUmRecords = "<<numRecords<<endl;
    int num_atoms = 0;
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();	
      ss << line;
      while(ss>>str) {
        if(str.size() > 4) {
          //multiple atom names together
          PSF::Atom temp_atom;
          unsigned int p=0, q=4;
          string sstr;
          while(p<str.size()) {
        temp_atom.number = num_atoms++;
        sstr = str.substr(p,q);
        temp_atom.atom_name = sstr;
        psf.atoms.push_back(temp_atom);
        p = p+q;
          }
        }
        else {
          PSF::Atom temp_atom;
          temp_atom.number = num_atoms++;
          temp_atom.atom_name = str;
          psf.atoms.push_back(temp_atom);
        }
                    
      }

      line = getline();
    }
    if(!(num_atoms ==  numRecords)) return false;
    //else cout<<"NUmATOMS = "<<num_atoms<<endl;
    //cout<<"PSF ATOMS ="<<psf.atoms.size()<<endl;
    //break;
      }
      else if (equalNocase("CHARGE",keyword)) {
    //line = getline();
    int index = 0;
        while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    std::vector<PSF::Atom>::iterator atom_iterator = psf.atoms.begin();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while(ss>>str) {
        (*atom_iterator).charge = toReal(str);
        index++;
        if(atom_iterator != psf.atoms.end()) atom_iterator++;
        else return false; //more entries for charge than # of atoms.
      }
      line = getline();
    }
    //cout<<"Index:Charge = "<<index<<endl;
    //break;

      }
      else if (equalStartNocase("MASS",keyword)) {
    //line = getline();
    int index = 0;
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    std::vector<PSF::Atom>::iterator atom_iterator = psf.atoms.begin();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while(ss>>str) {
        (*atom_iterator).mass = toReal(str);
        index++;
        if(atom_iterator != psf.atoms.end()) atom_iterator++;
        else return false; //more entries for mass than # of atoms.
      }
      line = getline();
    }
    //cout<<"Index:Mass = "<<index<<endl;
    //break;
      }
      else if(equalStartNocase("ATOM_TYPE_INDEX",keyword)) {
    //line = getline();
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1)!="%")) {
      ss.clear();
      ss << line;
      while(ss>>str) {
        iac.push_back(toInt(str));
      }
      line = getline();
    }
    //cout<<"IAC.size = "<<iac.size()<<endl;
    //break;
      }
      else if(equalStartNocase("NUMBER_EXCLUDED_ATOMS",keyword)) {
        //line = getline();
        do{
      line = getline();
        }while(line.substr(0,1) != "%");
      }
      else if(equalNocase("NONBONDED_PARM_INDEX",keyword)) {
    //line = getline();
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1)!="%")) {
      ss.clear();
      ss << line;
      while(ss>>str) {
        ico.push_back(toInt(str));
      }
      line = getline();
    }
    //break; 
      }
      else if(equalStartNocase("RESIDUE_LABEL",keyword)) {
    //int index = 0;
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    //std::vector<PSF::Atom>::iterator atom_iterator = psf.atoms.begin();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear(); 
      ss << line;
      while(ss>>str) {
        //(*atom_iterator).residue_name = str;
        //index++;
        //if(atom_iterator != psf.atoms.end()) atom_iterator++;
        //else return false; //more entries for charge than # of atoms.
        residues.push_back(str);	
      }
      line = getline();
    }
    //cout<<"Residues = "<<index<<endl;
    //cout<<"Total Residues = "<<residues.size()<<endl;
    //break;
      }
      else if(equalStartNocase("RESIDUE_POINTER",keyword)) {
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    std::vector<PSF::Atom>::iterator atom_iterator = psf.atoms.begin();
    std::vector<string>::iterator residue_label_iterator = residues.begin();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear(); 
      ss << line;
      while(ss>>str) {
        residue_pointers.push_back(toInt(str));
      }
      line = getline();		 
    }
    //std::vector<int>::iterator residue_ptr_iterator = residue_pointers.begin();
    int current, next, k = 0;

    while( (k+1) < ptrVals[11]){
      k++;
      current = residue_pointers[k-1];
      if ( (k-1) != (ptrVals[11]-1)) next = residue_pointers[k];
      else next = ptrVals[0];
      for(int p = current;p<next;p++){
        (*atom_iterator).residue_name = (*residue_label_iterator);
        atom_iterator++;
      }
      residue_label_iterator++;
    }
    //break;	
      }
      else if(equalStartNocase("BOND_FORCE_CONSTANT",keyword)) {
    //line = getline();
        while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while(ss >> str) {
        rk.push_back(toReal(str));
      }
      line = getline();
    }
    //cout<<"Bond Force Constants = "<<rk.size()<<endl;
    //break;
      }
      else if(equalStartNocase("BOND_EQUIL_VALUE",keyword)) {
    //line = getline();
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();	
      ss << line;
      while(ss >> str) {
        req.push_back(toReal(str));
      }
      line = getline();
        }
    //cout<<"Bond equil val = "<<req.size()<<endl;
    //break;	
      }
      else if(equalStartNocase("ANGLE_FORCE_CONSTANT",keyword)) {
    //line = getline();
    numAngles = ptrVals[4]+ptrVals[5];
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();	
      ss << line;
      while(ss>>str) {
        tk.push_back(toReal(str));
      }
      line = getline();
    }

    //cout<<"Angle Force Constant = "<<tk.size()<<endl;
    //break;	
      }
      else if(equalStartNocase("ANGLE_EQUIL_VALUE",keyword)) {
    //line = getline();
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while(ss >> str) {
        teq.push_back(toReal(str));
      }
      line = getline();
    }
    //cout<<"Angle Equil. Val = "<<teq.size()<<endl;
    //break;
      }
      else if(equalStartNocase("DIHEDRAL_FORCE_CONSTANT",keyword)) {
    //line = getline();
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while(ss >> str) {
        phi_k.push_back(toReal(str));
      }
      line = getline();
    }
    //cout<<"Dihedral_param_size = "<<phi_k.size()<<endl;
    //break;
      }
      else if(equalStartNocase("DIHEDRAL_PERIODICITY",keyword)) {
    //line = getline();
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();	
      ss << line;
      while(ss >> str) {
        phi_n.push_back(toReal(str));
      }
      line = getline();
    }
    //cout<<"Dihedral_param_size = "<<phi_n.size()<<endl;
    //break;
      }
      else if(equalStartNocase("DIHEDRAL_PHASE",keyword)) {
    //line = getline();
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();

    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while(ss >> str) {
        phi_phase.push_back(toReal(str));
      }
      line = getline();
    }
    //cout<<"Dihedral_param_size = "<<phi_phase.size()<<endl;
    //break;
      }
      else if(equalStartNocase("SOLTY",keyword)) {
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    do {
      line = getline();
    }while(line.substr(0,1) != "%");
    //cout<<line<<endl;
    //break;
      }
      else if(equalStartNocase("LENNARD_JONES_ACOEF",keyword)) {
                                        
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while( ss >> str ) {
        cn1.push_back(toReal(str));
      }
      line = getline();
      //cout<<line<<endl;
    }
    //cout<<"LJ param = "<<cn1.size()<<endl;
    //break;
      }
      else if(equalStartNocase("LENNARD_JONES_BCOEF",keyword)) {

    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while( ss >> str ) {
        cn2.push_back(toReal(str));
      }
      line = getline();
        }
    //cout<<"LJ param = "<<cn2.size()<<endl;
    //break;

      }
      else if(equalStartNocase("BONDS_INC_HYDROGEN",keyword)){
    //int temp_atom_number;
    //line = getline();
                                                       
    for(int i=0;i<ptrVals[2];i++) {
      int rkeq_index;
      PSF::Bond psf_bond;
      PAR::Bond par_bond;
                                                                                                            
      psf_bond.number = i+1;
      par_bond.number = i+1;
      file >> temp_atom_number;
      psf_bond.atom1 = (temp_atom_number+3)/3;
      file >> temp_atom_number;
      psf_bond.atom2 = (temp_atom_number+3)/3;
      file >> rkeq_index;

      par_bond.atom1 = psf_bond.atom1;
      par_bond.atom2 = psf_bond.atom2;
      par_bond.forceConstant = rk[rkeq_index];
      par_bond.distance = req[rkeq_index];
      psf.bonds.push_back(psf_bond);
      par.bonds.push_back(par_bond);
    }
    //cout<<"PSF Bonds = "<<psf.bonds.size()<<endl;
    //cout<<"PAR Bonds = "<<par.bonds.size()<<endl;
    line = getline();
    line = getline();
    //cout<<line<<endl;
    if(line.substr(0,1) != "%") {
      report <<recoverable
         <<"Mismatch in total number of bonds with hydrogens"
         <<endr;
      file.setstate(std::ios::failbit);
      file.close();
      return false;
    }
    //break;
      }
      else if(equalStartNocase("BONDS_WITHOUT_HYDROGEN",keyword)){
    //line = getline();

    for(int i=0;i<ptrVals[3];i++) {
      int rkeq_index;
      PSF::Bond psf_bond;
      PAR::Bond par_bond;

      psf_bond.number = i+1+ptrVals[15];
      par_bond.number = i+1+ptrVals[15];
      file >> temp_atom_number;
      psf_bond.atom1 = (temp_atom_number+3)/3;
      file >> temp_atom_number;
      psf_bond.atom2 = (temp_atom_number+3)/3;
      file >> rkeq_index;

      par_bond.atom1 = psf_bond.atom1;
      par_bond.atom2 = psf_bond.atom2;
      par_bond.forceConstant = rk[rkeq_index];
      par_bond.distance = req[rkeq_index];
      psf.bonds.push_back(psf_bond);
      par.bonds.push_back(par_bond);
    }

    //cout<<"PSF Bonds = "<<psf.bonds.size()<<endl;
    //cout<<"PAR Bonds = "<<par.bonds.size()<<endl;
    line = getline();
    line = getline();
    //cout<<line<<endl;
    if(line.substr(0,1) != "%") {
      report <<recoverable
         <<"Mismatch in total number of bonds without hydrogens"
         <<endr;
      file.setstate(std::ios::failbit);
      file.close();
      return false;
    }
      }
      else if(equalStartNocase("ANGLES_INC_HYDROGEN",keyword)) {

    for(int i=0;i<ptrVals[4];i++) {

      int tkeq_index;
      PSF::Angle psf_angle;
      PAR::Angle par_angle;

      psf_angle.number = i+1;
      par_angle.number = i+1;
      file >> temp_atom_number;
      psf_angle.atom1 = (temp_atom_number+3)/3;
      file >> temp_atom_number;
      psf_angle.atom2 = (temp_atom_number+3)/3;
      file >> temp_atom_number;
      psf_angle.atom3 = (temp_atom_number+3)/3;
      file >> tkeq_index;

      par_angle.atom1 = psf_angle.atom1;
      par_angle.atom2 = psf_angle.atom2;
      par_angle.atom3 = psf_angle.atom3;
      par_angle.forceConstant = tk[tkeq_index];
      par_angle.angleval = teq[tkeq_index];
      par_angle.ub_flag = 0;
      psf.angles.push_back(psf_angle);
      par.angles.push_back(par_angle);
    }
    line = getline();
    line = getline();
    if(line.substr(0,1) != "%") {
      report <<recoverable
         <<"Mismatch in total number of angles with hydrogens"
         <<endr;
      file.setstate(std::ios::failbit);
      file.close();
      return false;
    }

      }
      else if(equalStartNocase("ANGLES_WITHOUT_HYDROGEN",keyword)) {

    for(int i=0;i<ptrVals[5];i++) {

      int tkeq_index;
      PSF::Angle psf_angle;
      PAR::Angle par_angle;

      psf_angle.number = i+1+ptrVals[4];
      par_angle.number = i+1+ptrVals[5];
      file >> temp_atom_number;
      psf_angle.atom1 = (temp_atom_number+3)/3;
      file >> temp_atom_number;
      psf_angle.atom2 = (temp_atom_number+3)/3;
      file >> temp_atom_number;
      psf_angle.atom3 = (temp_atom_number+3)/3;
      file >> tkeq_index;

      par_angle.atom1 = psf_angle.atom1;
      par_angle.atom2 = psf_angle.atom2;
      par_angle.atom3 = psf_angle.atom3;
      par_angle.forceConstant = tk[tkeq_index];
      par_angle.angleval = teq[tkeq_index];
      par_angle.ub_flag = 0;
      psf.angles.push_back(psf_angle);
      par.angles.push_back(par_angle);
        }
        line = getline();
        line = getline();
    if(line.substr(0,1) != "%") {
      report <<recoverable
         <<"Mismatch in total number of angles without hydrogens"
         <<endr;
      file.setstate(std::ios::failbit);
      file.close();
      return false;
    }
    //break;
      }
      else if(equalStartNocase("DIHEDRALS_INC_HYDROGEN",keyword)) {

    for(int i=0;i<ptrVals[6];i++) {
      int phi_index;
      PSF::Dihedral psf_dihedral;
      PAR::Dihedral par_dihedral;

      psf_dihedral.number = i+1;
      par_dihedral.number = i+1;
      file >> temp_atom_number;
      psf_dihedral.atom1 = (temp_atom_number+3)/3;
      file >> temp_atom_number;
      psf_dihedral.atom2 = (temp_atom_number+3)/3;
      file >> temp_atom_number;
      psf_dihedral.atom3 = (temp_atom_number+3)/3;
      file >> temp_atom_number;
      psf_dihedral.atom4 = (temp_atom_number+3)/3;
      file >> phi_index;

      par_dihedral.atom1 = psf_dihedral.atom1;
      par_dihedral.atom2 = psf_dihedral.atom2;
      par_dihedral.atom3 = psf_dihedral.atom3;
      par_dihedral.atom4 = psf_dihedral.atom4;
      par_dihedral.multiplicity = 1;
      par_dihedral.forceConstant.push_back(phi_k[phi_index]);
      par_dihedral.periodicity.push_back((int)phi_n[phi_index]);
      par_dihedral.phaseShift.push_back(phi_phase[phi_index]);

      psf.dihedrals.push_back(psf_dihedral);
      par.dihedrals.push_back(par_dihedral);

    }
    line = getline();
    line = getline();
    if(line.substr(0,1) != "%") {
      report <<recoverable
         <<"Mismatch in total number of dihedrals with hydrogens"
         <<endr;
      file.setstate(std::ios::failbit);
      file.close();
      return false;
    }
      }
      else if(equalStartNocase("DIHEDRALS_WITHOUT_HYDROGEN",keyword)) {

        for(int i=0;i<ptrVals[7];i++) {
      int phi_index;
      PSF::Dihedral psf_dihedral;
      PAR::Dihedral par_dihedral;

      psf_dihedral.number = i+1+ptrVals[6];
      par_dihedral.number = i+1+ptrVals[6];
      file >> temp_atom_number;
      psf_dihedral.atom1 = abs((temp_atom_number+3)/3);
      file >> temp_atom_number;
      psf_dihedral.atom2 = abs((temp_atom_number+3)/3);
      file >> temp_atom_number;
      psf_dihedral.atom3 = abs((temp_atom_number+3)/3);
      file >> temp_atom_number;
      psf_dihedral.atom4 = abs((temp_atom_number+3)/3);
      file >> phi_index;

      par_dihedral.atom1 = psf_dihedral.atom1;
      par_dihedral.atom2 = psf_dihedral.atom2;
      par_dihedral.atom3 = psf_dihedral.atom3;
      par_dihedral.atom4 = psf_dihedral.atom4;

      par_dihedral.multiplicity = 0;
      phi_index--;
      do {
        par_dihedral.multiplicity++;
        phi_index++;
        par_dihedral.forceConstant.push_back(phi_k[phi_index]);
        par_dihedral.periodicity.push_back((int)phi_n[phi_index]);
        par_dihedral.phaseShift.push_back(phi_phase[phi_index]);
      }while(phi_n[phi_index]<0);

      psf.dihedrals.push_back(psf_dihedral);
      par.dihedrals.push_back(par_dihedral);
        }
    line = getline();
    line = getline();
    if(line.substr(0,1) != "%") {
      report <<recoverable
         <<"Mismatch in total number of dihedrals without hydrogens"
         <<endr;
      file.setstate(std::ios::failbit);
      file.close();
      return false;
    }
      }
      else if(equalStartNocase("EXCLUDED_ATOM_LIST",keyword)) {
    //Dont know what to do with these fields.
    while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    do{
      line = getline();
    }while(std::find(line.begin(), line.end(),'%') != line.begin());
      }
      else if(equalStartNocase("HBOND_ACOEF",keyword)) {
                                        
    //while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
    while((!line.empty()) && (std::find(line.begin(), line.end(),'%') == line.begin())){
      line = getline();
      //cout<<"In HBOND_ACOEF "<<line<<endl;
    }
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      //cout<<line<<endl;
      ss.clear();
      ss << line;
      while( ss >> str ) {
        hbond_acoff.push_back(toReal(str));
      }   
      line = getline();
      //cout<<line<<endl;
    }   
    //cout<<"HBOND_ACOEF = "<<hbond_acoff.size()<<endl;
    //break;
      }
      else if(equalStartNocase("HBOND_BCOEF",keyword)) {
                                        
    //while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
    while((!line.empty()) && (std::find(line.begin(), line.end(),'%') == line.begin()))
      line = getline();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while( ss >> str ) {
        hbond_bcoff.push_back(toReal(str));
      }
      line = getline();
    }

    //cout<<"HBOND_BCOEF = "<<hbond_bcoff.size()<<endl;
    //break;
      }
      else if(equalStartNocase("HBCUT", keyword)) {
    //ignoring this data
    //while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
    while((!line.empty()) && (std::find(line.begin(), line.end(),'%') == line.begin()))
      line = getline();
    do{
      line = getline();
    }while(std::find(line.begin(), line.end(),'%') != line.begin());

      }			  
      else if(equalStartNocase("AMBER_ATOM_TYPE", keyword)) {
    //fill up the atom types
    int index = 0;
        while((std::find(line.begin(), line.end(),'%') == line.begin()) || (line.empty()))
      line = getline();
    std::vector<PSF::Atom>::iterator atom_iterator = psf.atoms.begin();
    while((!line.empty()) && (line.substr(0,1) != "%")) {
      ss.clear();
      ss << line;
      while(ss>>str) {
        (*atom_iterator).atom_type = str;
        index++;
        if(atom_iterator != psf.atoms.end()) atom_iterator++;
        else return false; //more entries for charge than # of atoms.
      }
      line = getline();
    }
    cout<<"Total number of atoms = "<<index<<endl;
    break;
    //How to fill up PAR structure?

      } 

      else {
    line = getline();
      }
       


    }
    close();

    //Now open the parm99/gaff file to read in nonbonded parameters.
    if(!file.is_open()) {
        file.open("parm99.dat");
    }
    else {
        cout<<"Error opening parameter file"<<endl;
    }
    
    while(!file.eof()) {
        line = getline();
        ss << line;
        ss>>str;
        if(equalStartNocase("MOD4", str)) break;
    }
    
    int nonbonded_number = 0;
    while(!file.eof()) {
        PAR::Nonbonded nonb;
        line=getline();
        if(line.size() > 1) {
            ss << line;
            ss >> str;
            nonb.number = nonbonded_number;
            nonbonded_number++;
            nonb.atom = str;
            nonb.polarizability = 0;
            ss >> str;
            nonb.epsilon = toReal(str);
            ss >> str;
            nonb.sigma = toReal(str);
            nonb.negative = true;
            nonb.vdw = true;
            nonb.polarizability2 = 0.0;
            nonb.epsilon14 = nonb.epsilon;
            nonb.sigma14 = nonb.sigma;
            nonb.negative2 = true;
            //Insert into par structure
            par.nonbondeds.push_back(nonb);
        }
    }	
            
    return true;

  }
        
    

  PSF* AMBERReader::orphanPSF(){
    PSF* tmp = myPSF;
    myPSF = NULL;
    return tmp;
  }

  PAR* AMBERReader::orphanPAR(){
    PAR* tmp = myPAR;
    myPAR = NULL;
    return tmp;
  }

  void AMBERReader::writeData() {
    std::cout<<myPSF->atoms.size()<<std::endl;
    std::cout<<myPSF->bonds.size()<<std::endl;
    std::cout<<myPSF->angles.size()<<std::endl;
    std::cout<<myPAR->angles.size()<<std::endl;
    std::cout<<myPSF->dihedrals.size()<<std::endl;

  }

  AMBERReader& operator>>(AMBERReader& topReader, PSF& psf){
    if(topReader.myPSF == NULL){
      if(topReader.myPAR == NULL)
    topReader.myPAR = new PAR();
      topReader.read(psf, *(topReader.myPAR));
    }
    else {
      psf = *(topReader.myPSF);
      delete topReader.myPSF;
      topReader.myPSF = NULL;
    }
    return topReader;
  }

  AMBERReader& operator>>(AMBERReader& topReader, PAR& par){
    if(topReader.myPAR == NULL){
      if(topReader.myPSF == NULL)
    topReader.myPSF = new PSF();
      topReader.read(*(topReader.myPSF),par);
    }
    else {
      par = *(topReader.myPAR);
      delete topReader.myPAR;
      topReader.myPAR = NULL;
    }
    return topReader;
  }
}
