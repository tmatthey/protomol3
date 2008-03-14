#include <protomol/module/IOModule.h>

#include <protomol/io/PosVelReader.h>
#include <protomol/io/PSFReader.h>
#include <protomol/io/PARReader.h>

#include <protomol/config/Configuration.h>
#include <protomol/base/Exception.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/type/String.h>
#include <protomol/type/PDB.h>
#include <protomol/module/MainModule.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;

defineInputValue(InputVelocities, "velfile");
defineInputValueWithAliases(InputPositions, "posfile",
  ("coords")("coordinates"));
defineInputValueWithAliases(InputPSF, "psffile", ("structure"));
defineInputValueWithAliases(InputPAR, "parfile", ("parameters"));
defineInputValue(InputPDBScaling, "pdbScaling");
defineInputValue(InputDihedralMultPSF, "dihedralMultPSF");


void IOModule::init(ProtoMolApp *app) {
  Configuration *config = &app->config;
  
  InputPositions::registerConfiguration(config);
  InputVelocities::registerConfiguration(config);
  InputPSF::registerConfiguration(config);
  InputPAR::registerConfiguration(config);
  InputPDBScaling::registerConfiguration(config);
  InputDihedralMultPSF::registerConfiguration(config);
}

void IOModule::read(ProtoMolApp *app) {
  Configuration &config = app->config;

  // Positions
  PosVelReader reader;
  if (!reader.open(config[InputPositions::keyword]))
    THROW(string("Can't open position file '") +
      config[InputPositions::keyword].getString() + "'.");

  if (reader.tryFormat(PosVelReaderType::PDB)) {
    PDB pdb;
    if (!(reader >> pdb))
      THROW(string("Could not parse PDB position file '") +
        config[InputPositions::keyword].getString() + "'.");

    swap(app->positions, pdb.coords);

    // Add to output cache
    app->outputCache.add(pdb.atoms);

  } else if (!(reader >> app->positions))
    THROW(string("Could not parse position file '") +
      config[InputPositions::keyword].getString() +
      "'. Supported formats are : " +
      PosVelReaderType::getPossibleValues(", ") + ".");

  report << plain << "Using " << reader.getType() << " posfile '"
         << config[InputPositions::keyword] << "' ("
         << app->positions.size() << ")." << endr;

  // Velocities
  if (config.valid(InputVelocities::keyword)) {
    if (!reader.open(config[InputVelocities::keyword]))
      THROW(string("Can't open velocity file '") +
        config[InputVelocities::keyword].getString() + "'.");

    if (!(reader >> app->velocities))
      THROW(string("Could not parse velocity file '") +
        config[InputVelocities::keyword].getString() +
        "'. Supported formats are : " +
        PosVelReaderType::getPossibleValues(", ") + ".");

    report << plain << "Using " << reader.getType() << " velfile '"
           << config[InputVelocities::keyword] << "' ("
           << app->velocities.size() << ")." << endr;

    if (reader.getType() == "PDB" && (bool)config[InputPDBScaling::keyword]) {
      for (unsigned int i = 0; i < app->velocities.size(); i++)
        app->velocities[i] /= PDBVELSCALINGFACTOR;

      report << plain << "PDB velocities scaled." << endr;
    }

  } else if (config.valid(InputTemperature::keyword)) {
    app->velocities.resize(app->positions.size());

    report << plain << "Using temperature "
           << config[InputTemperature::keyword] << "K for the velocities  ("
           << app->velocities.size() << ")." << endr;
    // Create velocities later, we need the topology for that ...

  } else THROW("Neither temperature nor velocity file specified.");

  // PSF
  PSFReader psfReader;
  if (!psfReader.open(config[InputPSF::keyword]))
    THROW(string("Can't open PSF file '") +
      config[InputPSF::keyword].getString() + "'.");

  if (!(psfReader >> app->psf))
    THROW(string("Could not parse PSF file '") +
      config[InputPSF::keyword].getString() + "'.");

  report << plain << "Using PSF file '" << config[InputPSF::keyword]
         << "' (" << app->psf.atoms.size() << ")." << endr;

  // PAR
  PARReader parReader;
  if (!parReader.open(config[InputPAR::keyword]))
    THROW(string("Can't open PAR file '") +
      config[InputPAR::keyword].getString() + "'.");

  if (!(parReader >> app->par))
    THROW(string("Could not parse PAR file '") +
      config[InputPAR::keyword].getString() + "'.");

  report << plain << "Using PAR file '" << config[InputPAR::keyword]
         << "', " << (parReader.getCharmmTypeDetected() != PAR::CHARMM28 ?
                      "old" : "new") << " charmm force field.";

  if (!config[InputDihedralMultPSF::keyword].valid())
    config[InputDihedralMultPSF::keyword] =
      (parReader.getCharmmTypeDetected() != PAR::CHARMM28);

  if (config[InputDihedralMultPSF::keyword])
    report << " Dihedral multiplictity defined by PSF.";
  report << endr;

  // Test input
  if (app->positions.size() != app->velocities.size() ||
      app->positions.size() != app->psf.atoms.size())
    THROW("Positions, velocities and PSF input have different number "
          "of atoms.");  
}
