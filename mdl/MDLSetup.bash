#! /bin/bash

# USER: CHANGE THIS ACCORDINGLY
export PYTHONDIR=/afs/nd.edu/user25/tcickovs/Research/PYTHON2.5


export NUMPYDIR=${PYTHONDIR}
export GNUPLOTDIR=${PYTHONDIR}
export PYVERSION=python2.5
export TMPPATH=${PWD}:${PWD}/interface/:${PWD}/interface/integrator:${PWD}/interface/integrator/base/:${PWD}/interface/integrator/hessian/:${PWD}/interface/integrator/normal:${PWD}/interface/integrator/leapfrog:${PWD}/interface/topology:${PWD}/interface/io/:${PWD}/interface/output:${PWD}/interface/force:${PWD}/interface/force/system:${PWD}/interface/force/bonded:${PWD}/interface/force/nonbonded:${PWD}/interface/type:${PWD}/interface/base:${PWD}/src/:${PWD}/src/factories:${PWD}/src/modifiers:${PWD}/src/propagators:${PWD}/src/propagators/objects:${PWD}/src/toplevel:${GNUPLOTDIR}/lib/${PYVERSION}/site-packages/Gnuplot:${NUMPYDIR}/lib/${PYVERSION}/site-packages/:${NUMPYDIR}/lib/${PYVERSION}/site-packages/numpy/

if [[ "$PYTHONPATH" != "" ]] 
then export PYTHONPATH=${PYTHONDIR}/lib/${PYVERSION}:${TMPPATH}:${PYTHONPATH}
else export PYTHONPATH=${PYTHONDIR}/lib/${PYVERSION}:${TMPPATH}
fi 
 
export MDLROOT=${PWD}
export PATH=${PYTHONDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${PWD}/interface/:${PYTHONDIR}/lib/${PYVERSION}:${LD_LIBRARY_PATH}
export PROTOMOL_HOME=${PWD}/../

