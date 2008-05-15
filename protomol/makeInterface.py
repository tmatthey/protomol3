# These modules should not have interface files generated
# i.e, there are special 'rules' applied (for example for template
# instantiatations.
excluded_modules = ['GenericTopology', 'PARReader', 'PSFReader', 'PDBReader', 'XYZReader', 'XYZBinReader', 'EigenvectorReader', 'EigenvectorTextReader', 'Vector3DBlock', 'TopologyUtilities', 'ForceGroup', 'OutputCache', 'PDBWriter', 'XYZWriter', 'ScalarStructure']

def exclude(module):
    return excluded_modules.count(module)

def makeInterface(dir, module):
    filename = dir+'/'+module+'.i'
    FILE = open(filename,"w")
    FILE.writelines('%module '+module+'\n')
    FILE.writelines('%{\n')
    FILE.writelines('#include \"'+module+'.h\"\n')
    FILE.writelines('#include <protomol/type/Real.h>\n')
    FILE.writelines('#include <protomol/ProtoMolApp.h>\n')
    if (module == 'ProtoMolApp'):
       FILE.writelines('#include <protomol/integrator/STSIntegrator.h>\n')
    if (dir.find('integrator/') != -1):
       FILE.writelines('#include <protomol/integrator/Integrator.h>\n')
       FILE.writelines('#include <protomol/integrator/StandardIntegrator.h>\n')
       FILE.writelines('#include <protomol/integrator/STSIntegrator.h>\n')
       FILE.writelines('#include <protomol/integrator/MTSIntegrator.h>\n')
       if (dir.find('normal') != -1):
           FILE.writelines('#include <protomol/type/EigenvectorInfo.h>\n')
       FILE.writelines('#include \"ndarrayobject.h\"\n')
    elif (dir.find('output') != -1 and module != 'OutputCache'):
       FILE.writelines('#include <protomol/output/Output.h>\n')
       FILE.writelines('#include <protomol/output/OutputFile.h>\n')
       FILE.writelines('#include <protomol/output/OutputCache.h>\n')
    FILE.writelines('using namespace ProtoMol;\n')
    FILE.writelines('%}\n')
    FILE.writelines('\n')
    FILE.writelines('%include <protomol/type/Real.h>\n')
    if (dir.find('integrator/') != -1):
       FILE.writelines('%feature (\"dynamic_cast\");\n')
       FILE.writelines('%include \"std_string.i\"\n')
       FILE.writelines('%include <protomol/integrator/Integrator.h>\n')
       FILE.writelines('%include <protomol/integrator/StandardIntegrator.h>\n')
       FILE.writelines('%include <protomol/integrator/STSIntegrator.h>\n')
       FILE.writelines('%include <protomol/integrator/MTSIntegrator.h>\n')
       if (dir.find('normal') != -1):
           FILE.writelines('%include <protomol/type/EigenvectorInfo.h>\n')
           FILE.writelines('%include <protomol/integrator/normal/NormalModeUtilities.h>\n')
    elif (dir.find('output') != -1 and module != 'OutputCache'):
       FILE.writelines('\n')
       FILE.writelines('%include \"std_string.i\"\n')
       FILE.writelines('%include <protomol/output/Output.h>\n')
       FILE.writelines('%include <protomol/output/OutputFile.h>\n')
       FILE.writelines('%include <protomol/output/OutputCache.h>\n')
    elif (module.find('Writer') != -1):
       FILE.writelines('%include \"std_string.i\"\n')
       FILE.writelines('%include \"std_vector.i\"\n')
       FILE.writelines('%template() std::vector<ProtoMol::Atom>;\n')
       FILE.writelines('%template() std::vector<ProtoMol::AtomType>;\n')
       FILE.writelines('%include <protomol/type/Vector3DBlock.i>\n')
       FILE.writelines('%include <protomol/topology/Atom.h>\n')
       FILE.writelines('%include <protomol/topology/AtomType.h>\n')
       FILE.writelines('%include <protomol/io/File.h>\n')
       FILE.writelines('%include <protomol/io/Writer.h>\n')
    elif (module.find('Reader') != -1):
       FILE.writelines('%include \"std_string.i\"\n')
    elif (module == 'ProtoMolApp'):
       FILE.writelines('%include <protomol/type/Vector3DBlock.i>\n')
       FILE.writelines('%include <protomol/type/ScalarStructure.h>\n')
    FILE.writelines('%include \"'+module+'.h\"\n')
    if (dir.find('integrator/') != -1 and module != 'NormalModeUtilities'):
       FILE.writelines('\n')
       FILE.writelines('%extend ProtoMol::'+module+' {\n')
       FILE.writelines('ProtoMolApp* appInit(GenericTopology* topo,\n')
       FILE.writelines('                     Vector3DBlock& positions,\n')
       FILE.writelines('                      Vector3DBlock& velocities,\n')
       FILE.writelines('                      ScalarStructure energies) {\n')
       FILE.writelines('   import_array1(NULL);\n')
       FILE.writelines('   ProtoMolApp* app = new ProtoMolApp();\n')
       FILE.writelines('   app->topology = topo;\n')
       FILE.writelines('   app->positions.vec = positions.vec;\n')
       FILE.writelines('   app->positions.c = positions.c;\n')
       FILE.writelines('   app->velocities.vec = velocities.vec;\n')
       FILE.writelines('   app->velocities.c = velocities.c;\n')
       FILE.writelines('   app->energies = energies;\n')
       FILE.writelines('   self->initialize(app);\n')
       FILE.writelines('   app->integrator = self;\n')
       FILE.writelines('   app->outputCache.initialize(app);\n')
       FILE.writelines('   return app;\n')
       FILE.writelines('}\n')
       FILE.writelines('};\n')
    elif (dir.find('output') != -1 and module != 'OutputCache'):
       FILE.writelines('\n')
       FILE.writelines('%extend ProtoMol::'+module+' {\n')
       FILE.writelines('void uncache(ProtoMolApp* app) {\n')
       FILE.writelines('   app->outputCache.uncache();\n')
       FILE.writelines('}\n')
       FILE.writelines('};\n')
    elif (module == 'ProtoMolApp'):
       FILE.writelines('\n')
       FILE.writelines('%extend ProtoMol::'+module+' {\n')
       FILE.writelines('void makeApp(GenericTopology* topo,\n')
       FILE.writelines('             ProtoMol::Vector3DBlock positions,\n')
       FILE.writelines('             ProtoMol::Vector3DBlock velocities,\n')
       FILE.writelines('             ScalarStructure energies,\n')
       FILE.writelines('             Real timestep) {\n')
       FILE.writelines('   self->topology = topo;\n')
       FILE.writelines('   self->positions.vec = positions.vec;\n')
       FILE.writelines('   self->positions.c = positions.c;\n')
       FILE.writelines('   self->velocities.vec = velocities.vec;\n')
       FILE.writelines('   self->velocities.c = velocities.c;\n')
       FILE.writelines('   self->energies = energies;\n')
       #FILE.writelines('   self->initialize(app);\n')
       FILE.writelines('   self->integrator = new STSIntegrator(timestep, NULL);\n')
       FILE.writelines('   self->outputCache.initialize(self);\n')
       FILE.writelines('}\n')
       FILE.writelines('};\n')       

       
