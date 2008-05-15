#! /bin/bash

# USER: CHANGE THESE ACCORDINGLY
export PYTHONDIR=/afs/nd.edu/user25/tcickovs/Research/PYTHON2.5
export LAPACKDIR=/afs/nd.edu/user25/tcickovs/Research/LAPACK
export NUMPYDIR=/afs/nd.edu/user25/tcickovs/Research/PYTHON2.5
export GNUPLOTDIR=/afs/nd.edu/user25/tcickovs/Research/PYTHON2.5
export PYVERSION=python2.5

export TMPPATH=${PWD}:${PWD}/../protomol/integrator/base/:${PWD}/../protomol/integrator/hessian/:${PWD}/../protomol/integrator/normal:${PWD}/../protomol/integrator/leapfrog:${PWD}/../protomol/topology:${PWD}/../protomol/io/:${PWD}/../protomol/output:${PWD}/../protomol/force:${PWD}/../protomol/force/system:${PWD}/../protomol/force/bonded:${PWD}/../protomol/force/nonbonded:${PWD}/../protomol/type:${PWD}/../protomol/base:${PWD}/src:${PWD}/src/factories:${PWD}/src/modifiers:${PWD}/src/propagators:${PWD}/src/propagators/objects:${PWD}/src/toplevel:${GNUPLOTDIR}/lib/${PYVERSION}/site-packages/Gnuplot:${NUMPYDIR}/lib/${PYVERSION}/site-packages/:${NUMPYDIR}/lib/${PYVERSION}/site-packages/numpy/

if [[ "$PYTHONPATH" != "" ]] 
then export PYTHONPATH=${TMPPATH}:${PYTHONPATH}
else export PYTHONPATH=${TMPPATH}
fi 
 
export MDLROOT=${PWD}
export PATH=${PYTHONDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${PYTHONDIR}/lib/${PYVERSION}:${LD_LIBRARY_PATH}
export PYTHONPATH=${PYTHONDIR}/lib/${PYVERSION}:${PYTHONPATH}
export LD_LIBRARY_PATH=${PWD}/../:${LAPACKDIR}/lib:${LD_LIBRARY_PATH}
export PYTHONPATH=${LAPACKDIR}/lib:${PYTHONPATH}

