#include <protomol/topology/CoulombSCPISMParameters.h>
#include <protomol/topology/CoulombSCPISMParameterTable.h>
#include <protomol/base/Report.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

void CoulombSCPISMParameterTable::populateTable() {
  CoulombSCPISMParameters temp;
  // Order:
  // alpha_i, hbond_factor, R_iw, sqrt_alpha_i, r_cov, gamma_i, isHbonded, A_i,
  // B_i, C_i, R_vdw
  
  temp.set(0.4906, 0.6259, 0.5000, 0.7004, 0.3700, 0.0052, PH, 1.20,
           3.68, 0.80, 0.224500);
  myData["H"] = temp; // polar H
  temp.set(0.5300, 9.6746, 0.5000, 0.7280, 0.3700, 0.0052, PH, 1.20,
           3.68, 0.80, 0.224500);
  myData["HC"] = temp; // N-ter H
  temp.set(0.5300, 0.5000, 0.5000, 0.7280, 0.3700, 0.0052, NO, 17.0, 
           9.00, 0.50, 1.320000);
  myData["HA"] = temp; // nonpolar H
  temp.set(0.4906, 2.3800, 0.5000, 0.7004, 0.3700, 0.0052, PH, 1.20, 
           3.68, 0.80, 0.224500);
  myData["HT"] = temp; // TIPS3P Water H
  temp.set(0.5274, 0.5000, 0.5000, 0.7262, 0.3700, 0.0052, NO, 17.0,  
           9.00, 0.50, 1.358200);
  myData["HP"] = temp; // aromatic H
  temp.set(0.4858, 0.5000, 0.5000, 0.6970, 0.3700, 0.0052, NO, 17.0,  
           9.00, 0.50, 1.320000);
  myData["HB"] = temp; // backbone H
  temp.set(0.4651, 0.5000, 0.5000, 0.6820, 0.3700, 0.0052, NO, 17.0,  
           9.00, 0.50, 0.900000);
  myData["HR1"] = temp; //his he1, (+) his HG,HD2
  temp.set(0.4580, 0.5000, 0.5000, 0.6768, 0.3700, 0.0052, NO, 17.0,  
           9.00, 0.50, 0.700000);
  myData["HR2"] = temp; //(+) his HE1
  temp.set(0.4893, 0.5000, 0.5000, 0.6995, 0.3700, 0.0052, NO, 17.0,  
           9.00, 0.50, 1.468000);
  myData["HR3"] = temp; //neutral his HG, HD2
  temp.set(0.4859, 0.5000, 0.5000, 0.6971, 0.3700, 0.0052, NO, 17.0,  
           9.00, 0.50, 0.450000);
  myData["HS"] = temp; // thiol H
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.3700, 0.0052, NO, 17.0,  
           9.00, 0.50, 0.224500);
  myData["HA1"] = temp; // for alkene; RHC=CR, wildcard Rvdw
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.3700, 0.0052, NO, 17.0,  
           9.00, 0.50, 0.224500);
  myData["HA2"] = temp; // for alkene; H2C=CR, wildcard Rvdw
  temp.set(0.5298, 0.5000, 0.5000, 0.7273, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.000000);
  myData["C"] = temp; // polar C
  temp.set(0.5192, 0.5000, 0.5000, 0.7206, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 1.992400);
  myData["CA"] = temp; // aromatic C
  temp.set(0.5209, 0.5000, 0.5000, 0.7217, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.275000);
  myData["CT1"] = temp; //aliphatic sp3 C for CH
  temp.set(0.5298, 0.5000, 0.5000, 0.7273, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.175000);
  myData["CT2"] = temp; //aliphatic sp3 C for CH2
  temp.set(0.5082, 0.5000, 0.5000, 0.7129, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.060000);
  myData["CT3"] = temp; //aliphatic sp3 C for CH3
  temp.set(0.4538, 0.5000, 0.5000, 0.6736, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 1.800000);
  myData["CPH1"] = temp; //his CG and CD2 carbons
  temp.set(0.4974, 0.5000, 0.5000, 0.7053, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 1.800000);
  myData["CPH2"] = temp; //his CE1 carbon
  temp.set(0.5217, 0.5000, 0.5000, 0.7229, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 1.800000);
  myData["CPT"] = temp; //trp C between rings
  temp.set(0.5240, 0.5000, 0.5000, 0.7239, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 1.992400);
  myData["CY"] = temp; //TRP C in pyrrole ring
  temp.set(0.4763, 0.5000, 0.5000, 0.6901, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.275000);
  myData["CP1"] = temp; //tetrahedral C (proline CA)
  temp.set(0.4599, 0.5000, 0.5000, 0.6782, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.175000);
  myData["CP2"] = temp; //tetrahedral C (proline CB/CG)
  temp.set(0.4627, 0.5000, 0.5000, 0.6802, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.175000);
  myData["CP3"] = temp; //tetrahedral C (proline CD)
  temp.set(0.4970, 0.5000, 0.5000, 0.7050, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.000000);
  myData["CC"] = temp; //carbonyl C for asn,asp,gln,glu
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.000000);
  myData["CD"] = temp; //carbonyl C for none amides, asp,glu,cter
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 1.800000);
  myData["CPA"] = temp; //heme alpha-C
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 1.800000);
  myData["CPB"] = temp; //heme beta-C
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 1.800000);
  myData["CPM"] = temp; //heme meso-C
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.100000);
  myData["CM"] = temp; //heme CO carbon
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.200000);
  myData["CS"] = temp; //thiolate carbon
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.090000);
  myData["CE1"] = temp; //for alkene; RHC=CR
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7700, 0.0052, NO, 31.0,  
           7.92, 0.32, 2.080000);
  myData["CE2"] = temp; //for alkene; H2C=CR
  temp.set(0.4894, 0.5000, 0.5000, 0.6996, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["N"] = temp; // proline N
  temp.set(0.4598, 0.5000, 0.5000, 0.6781, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NR1"] = temp; // neutral his protonated ring N
  temp.set(0.4839, 19.4824, 0.6956, 0.5000, 0.7400, 0.0052, PA, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NR2"] = temp; // neutral his unprotonated ring N
  temp.set(0.4742, 0.5000, 0.5000, 0.6886, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NR3"] = temp; // charged his  ring N
  temp.set(0.5103, 0.5000, 0.5000, 0.7144, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NH1"] = temp; // peptide N
  temp.set(0.5045, 0.5000, 0.5000, 0.7103, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NH2"] = temp; // amide N
  temp.set(0.5300, 0.5000, 0.5000, 0.7280, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NH3"] = temp; // ammonium N
  temp.set(0.5299, 0.5000, 0.5000, 0.7279, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NC2"] = temp; // guanidinium N
  temp.set(0.5284, 0.5000, 0.5000, 0.7269, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NY"] = temp; // TRP N in pyrrole ring
  temp.set(0.4894, 0.5000, 0.5000, 0.6996, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NP"] = temp; // Proline ring NH2+ (N-terminal)
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO, 60.0,  
           9.66, 0.22, 1.850000);
  myData["NPH"] = temp; //heme pyrrole N
  temp.set(0.4810, 4.0085, 0.5000, 0.6935, 0.7400, 0.0052, PA, 40.0,  
           8.52, 0.29, 1.700000);
  myData["O"] = temp; // carbonyl oxygen
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO, 40.0,  
           8.52, 0.29, 1.700000);
  myData["OB"] = temp; //carbonyl oxygen in acetic acid
  temp.set(0.6000, 8.1838, 0.5000, 0.7746, 0.7400, 0.0052, PA, 40.0,  
           8.52, 0.29, 1.700000);
  myData["OC"] = temp; // carboxlate oxygen
  temp.set(0.4843, 3.3433, 0.5000, 0.6959, 0.7400, 0.0052, PA, 40.0,  
           8.52, 0.29, 1.770000);
  myData["OH1"] = temp; // hydroxyl oxygen
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO, 40.0,  
           8.52, 0.29, 1.770000);
  myData["OS"] = temp; // ester oxygen
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO, 40.0,  
           8.52, 0.29, 1.768200);
  myData["OT"] = temp; // TIPS3P Water O
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.7400, 0.0052, NO, 40.0,  
           8.52, 0.29, 1.700000);
  myData["OM"] = temp; // heme CO/O2 oxygen
  temp.set(0.5054, 5.0037, 0.5000, 0.7109, 1.0400, 0.0052, PA, 60.0,  
           9.10, 0.22, 2.000000);
  myData["S"] = temp; // sulphur
  temp.set(0.5054, 0.5000, 0.5000, 0.7109, 1.0400, 0.0052, NO, 60.0,  
           9.10, 0.22, 1.975000);
  myData["SM"] = temp; // sulphur C-S-S-C type
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 1.0400, 0.0052, NO, 60.0,  
           9.10, 0.22, 2.200000);
  myData["SS"] = temp; // thiolate sulphur
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO, 0.0,  
           0.00, 0.00, 1.4800);
  myData["HE"] = temp; // helium
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO, 0.0,  
           0.00, 0.00, 1.5300);
  myData["NE"] = temp; // neon
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO, 0.0,  
           0.00, 0.00, 1.3670);
  myData["CAL"] = temp; // calcium 2+
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO, 0.0,  
           0.00, 0.00, 1.0900);
  myData["ZN"] = temp; // zinc (II) cation
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO, 0.0,  
           0.00, 0.00, 0.6500);
  myData["FE"] = temp; // heme iron 56
  temp.set(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, NO, 0.0,  
           0.00, 0.00, 0.0000);
  myData["DUM"] = temp; // dummy atom
}

void CoulombSCPISMParameterTable::displayTable() {
  report << "Coulomb SCPISM Parameter Table!!!" << endl;
  
  map<string, CoulombSCPISMParameters>::const_iterator it;
  for (it = myData.begin(); it != myData.end(); it++) {
    report << it->first << "\t";
    report << it->second.alpha_i << "\t";
    report << it->second.hbond_factor << "\t";
    report << it->second.R_iw << "\t";
    //report << it->second.sqrt_alpha_i << "\t";
    report << it->second.r_cov << "\t";
    report << it->second.gamma_i << "\t";
    report << it->second.isHbonded << "\t";
    report << it->second.A_i << "\t";
    report << it->second.B_i << "\t";
    report << it->second.C_i << "\t";
    //report << it->second.R_vdw << "\t";
    report << endl;
  }
}

