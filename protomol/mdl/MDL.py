import os
import sys

pythondir = sys.exec_prefix
numpydir = pythondir
gnuplotdir = pythondir
pyversion = 'python2.5'

pwd = os.getenv('PWD')

sys.path.append(pwd)
sys.path.append(pwd+'/interface/')
sys.path.append(pwd+'/interface/integrator')
sys.path.append(pwd+'/interface/integrator/base/')
sys.path.append(pwd+'/interface/integrator/hessian/')
sys.path.append(pwd+'/interface/integrator/normal/')
sys.path.append(pwd+'/interface/integrator/leapfrog/')
sys.path.append(pwd+'/interface/topology/')
sys.path.append(pwd+'/interface/io/')
sys.path.append(pwd+'/interface/output/')
sys.path.append(pwd+'/interface/force/')
sys.path.append(pwd+'/interface/force/system/')
sys.path.append(pwd+'/interface/force/bonded/')
sys.path.append(pwd+'/interface/force/nonbonded/')
sys.path.append(pwd+'/interface/type/')
sys.path.append(pwd+'/interface/base/')
sys.path.append(pwd+'/src/')
sys.path.append(pwd+'/src/factories/')
sys.path.append(pwd+'/src/modifiers/')
sys.path.append(pwd+'/src/propagators/')
sys.path.append(pwd+'/src/propagators/objects/')
sys.path.append(pwd+'/src/toplevel')
sys.path.append(gnuplotdir+'/lib/'+pyversion+'/site-packages/Gnuplot/')
sys.path.append(numpydir+'/lib/'+pyversion+'/site-packages/')
sys.path.append(numpydir+'/lib/'+pyversion+'/site-packages/numpy/')


os.environ['MDLROOT'] = pwd
os.system('export LD_LIBRARY_PATH='+pythondir+'/lib/'+pyversion+':$LD_LIBRARY_PATH')

if (not os.environ.has_key('PROTOMOL_HOME')):
    os.environ['PROTOMOL_HOME'] = pwd+':/../'


class Prompt:
  def __init__(self, str='>>>'):
    self.str = str;
        
  def __str__(self):
    return self.str
    
sys.ps1 = Prompt("\033[1;31m[MDL] \033[0m")


from Physical import *
from IO import *
from ForceField import *
from Forces import *
from Propagator import *


