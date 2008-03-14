#include <protomol/module/NormalModeModule.h>

#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/type/String.h>
#include <protomol/io/EigenvectorReader.h>
#include <protomol/io/EigenvectorTextReader.h>

#include <protomol/integrator/normal/NormalModeLangevin.h>
#include <protomol/integrator/normal/NormalModeMinimizer.h>
#include <protomol/integrator/normal/NormalModeDiagonalize.h>
#include <protomol/integrator/normal/NormalModeMori.h>
#include <protomol/integrator/normal/NormalModeRelax.h>
#include <protomol/integrator/normal/NormalModeBrownian.h>

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

defineInputValue(InputEigenVectors, "eigfile");
defineInputValue(InputEigTextFile, "eigtextfile");

void NormalModeModule::init(ProtoMolApp *app) {
  InputEigenVectors::registerConfiguration(&app->config);
  InputEigTextFile::registerConfiguration(&app->config);

  app->integratorFactory.registerExemplar(new NormalModeLangevin());
  app->integratorFactory.registerExemplar(new NormalModeMinimizer());
  app->integratorFactory.registerExemplar(new NormalModeDiagonalize());
  app->integratorFactory.registerExemplar(new NormalModeMori());
  app->integratorFactory.registerExemplar(new NormalModeRelax());
  app->integratorFactory.registerExemplar(new NormalModeBrownian());
}

void NormalModeModule::read(ProtoMolApp *app) {
  Configuration &config = app->config;

  // Eigenvectors/values
  if (config.valid(InputEigTextFile::keyword)) {
    EigenvectorTextReader evTextReader;

    if (config.valid(InputEigTextFile::keyword)) {
      if (!evTextReader.open(config[InputEigTextFile::keyword]))
        THROWS("Can't open eigenvector file '"
               << config[InputEigTextFile::keyword].getString() << "'.");

      if (!(evTextReader >> app->eigenInfo)) {
        THROWS("Could not parse eigenvector file '"
               << (string)config[InputEigTextFile::keyword] << "'");

        if (app->eigenInfo.myEigenvectorLength != (double)app->positions.size())
          THROWS("Eigenvector length is wrong, should be "
                 << app->positions.size() << " got "
                 << app->eigenInfo.myEigenvectorLength << ".");

        if (app->eigenInfo.myNumEigenvectors < 1 ||
            app->eigenInfo.myNumEigenvectors > (double)app->positions.size())
          THROWS("Wrong number of eigenvectors ("
                 << app->eigenInfo.myNumEigenvectors << ").");
      }

      report << plain << "Using eigfile '" << config[InputEigTextFile::keyword]
             << "' (" << app->eigenInfo.myEigenvectorLength << ")." << endr;
      eiValid = true;
    }

  } else if (config.valid(InputEigenVectors::keyword)) {
    EigenvectorReader evReader;

    if (!evReader.open(config[InputEigenVectors::keyword]))
      THROWS("Can't open eigenvector file '"
             << (string)config[InputEigenVectors::keyword] << "'.");

    if (!(evReader >> app->eigenInfo)) {
      if (app->eigenInfo.myEigenvectorLength != (double)app->positions.size())
        THROWS("Eigenvector length is wrong, should be "
               << app->positions.size() << " got "
               << app->eigenInfo.myEigenvectorLength << ".");

      if (app->eigenInfo.myNumEigenvectors < 1 ||
          app->eigenInfo.myNumEigenvectors > (double)app->positions.size())
        THROWS("Wrong number of eigenvectors ("
               << app->eigenInfo.myNumEigenvectors << ").");
    }

    report << plain << "Using eigfile '"
           << config[InputEigenVectors::keyword] << "' ("
           << app->eigenInfo.myEigenvectorLength << ")." << endr;
    
    eiValid = true;
  }
}

void NormalModeModule::postBuild(ProtoMolApp *app) {
  // Normal mode?
  // New method tries a dynamic cast, then updates the pointers from ei
  NormalModeUtilities *nmint; // dynamic cast working
  nmint  = dynamic_cast<NormalModeUtilities *>(app->integrator);
  if (nmint) { // cast worked?
	nmint->setIntegratorSetPointers(app->integrator, &app->eigenInfo, eiValid);
	report << plain << "Using new Normal Mode integrator. " << endr;

  } else if (eiValid)
    report << plain << "Warning: Eigenvector file defined but using "
           << "non-Normal Mode integrator!" << endr;
}
