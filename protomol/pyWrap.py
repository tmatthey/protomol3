import sys
import os
from makeInterface import *


# List of Python modules which should be generated.
# The mapping is from a directory in the ProtoMol framework
# to a list of module names.
# The module names will correspond to prefixes for the interface
# file, source file and shared object.
mdlmodules = {'protomol/integrator/leapfrog':['LeapfrogIntegrator', 'LeapfrogTruncatedShadow', 'PLeapfrogIntegrator', 'DMDLeapfrogIntegrator', 'NoseNVTLeapfrogIntegrator'],
              'protomol/integrator/':['STSIntegrator', 'MTSIntegrator'],
	      'protomol/integrator/base':['LangevinImpulseIntegrator', 'CGMinimizerIntegrator', 'NumericalDifferentiation'],
              'protomol/type':['Vector3DBlock', 'ScalarStructure'],
	      'protomol/integrator/normal':['NormalModeBrownian', 'NormalModeDiagonalize', 'NormalModeMinimizer', 'NormalModeLangevin', 'NormalModeUtilities', 'NormalModeMori', 'NormalModeRelax'],
              'protomol/integrator/hessian':['HessianInt'],
              'protomol/io':['DCDTrajectoryReader', 'EigenvectorReader', 'EigenvectorTextReader', 'PARReader', 'PDBReader', 'PDBWriter', 'PSFReader', 'XYZBinReader', 'XYZReader', 'XYZTrajectoryReader', 'XYZTrajectoryWriter', 'XYZWriter'],
              'protomol/output':['OutputCache', 'OutputDCDTrajectory', 'OutputDCDTrajectoryVel', 'OutputEnergies', 'OutputFinalPDBPos', 'OutputFinalXYZPos', 'OutputFinalXYZVel', 'OutputScreen', 'OutputXYZTrajectoryForce', 'OutputXYZTrajectoryPos', 'OutputXYZTrajectoryVel'],
              'protomol/base':['MathUtilities'],
	      'protomol/topology':['TopologyUtilities', 'GenericTopology'],
	      'protomol/force/bonded':['BondForce','AngleForce','DihedralForce','HarmDihedralForce','ImproperForce'],
              'protomol/force/nonbonded':['SimpleFullForce', 'CutoffForce'],
              'protomol/force/system':['PySystemForce'],
              'protomol/force':['ForceGroup'],
              'protomol':['ProtoMolApp']
	     }


def pyWrap(env):
  import distutils.sysconfig
  env.Append(SWIGFLAGS=['-c++', '-python', '-w312', '-w314', '-w315', '-w317', '-w361', '-w362', '-w389', '-w401', '-w454', '-w503', '-w509'],
                     CPPPATH=[distutils.sysconfig.get_python_inc()],
                     SHLIBPREFIX="",
                     ENV={'PATH':os.environ['PATH']})

  pyvers = "python%d.%d" % sys.version_info[:2]
  numpypath = sys.exec_prefix+"/lib/"+pyvers+"/site-packages/numpy/core/include/numpy/"

  #env['CXX'] = "icpc"
  env['ARFLAGS'] = "rcS"
  env.Append(CPPPATH = '#')
  env.Append(LINKFLAGS=' -Wl,-E')
  env.Append(CXXFLAGS=' -fPIC')
  env.Append(SHCXXFLAGS=' -I'+numpypath)
  env.Append(SHLINKFLAGS=' -Wl,-E')
  env.Append(_LIBDIRFLAGS="-L.")

  for dir in mdlmodules.iterkeys():
    modulelist = mdlmodules[dir]
    for i in range(0, len(modulelist)):
       module = modulelist[i]
       if (dir.find('force/') == -1 and excluded_modules.count(module) == 0):
          makeInterface(dir, module)
       env.SharedLibrary(target=dir+'/_'+module+'.so', source=[dir+'/'+module+'.i'], SHLIBPREFIX="")
