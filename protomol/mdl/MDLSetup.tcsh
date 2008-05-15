#! /bin/tcsh

# USER: CHANGE THESE ACCORDINGLY
setenv PYTHONDIR /afs/nd.edu/user25/tcickovs/Research/PYTHON2.5
setenv LAPACKDIR /afs/nd.edu/user25/tcickovs/Research/LAPACK
setenv NUMPYDIR /afs/nd.edu/user25/tcickovs/Research/PYTHON2.5
setenv GNUPLOTDIR /afs/nd.edu/user25/tcickovs/Research/PYTHON2.5
setenv PYVERSION python2.5

setenv TMPPATH ${PWD}:${PWD}/../protomol/:${PWD}/../protomol/integrator:${PWD}/../protomol/integrator/base/:${PWD}/../protomol/integrator/hessian/:${PWD}/../protomol/integrator/normal:${PWD}/../protomol/integrator/leapfrog:${PWD}/../protomol/topology:${PWD}/../protomol/io/:${PWD}/../protomol/output:${PWD}/../protomol/force:${PWD}/../protomol/force/system:${PWD}/../protomol/force/bonded:${PWD}/../protomol/force/nonbonded:${PWD}/../protomol/type:${PWD}/../protomol/base:${PWD}/src/:${PWD}/src/factories:${PWD}/src/modifiers:${PWD}/src/propagators:${PWD}/src/propagators/objects:${PWD}/src/toplevel:${GNUPLOTDIR}/lib/${PYVERSION}/site-packages/Gnuplot:${NUMPYDIR}/lib/${PYVERSION}/site-packages/:${NUMPYDIR}/lib/${PYVERSION}/site-packages/numpy/

if ($?PYTHONPATH != 0) then
   setenv PYTHONPATH ${TMPPATH}:${PYTHONPATH}
else
   setenv PYTHONPATH ${TMPPATH}
endif

setenv MDLROOT ${PWD}
setenv PATH ${PYTHONDIR}/bin:${PATH}
setenv LD_LIBRARY_PATH ${PYTHONDIR}/lib/${PYVERSION}:${LD_LIBRARY_PATH}
setenv PYTHONPATH ${PYTHONDIR}/lib/${PYVERSION}:${PYTHONPATH}
setenv LD_LIBRARY_PATH ${PWD}/../:${LAPACKDIR}/lib:${LD_LIBRARY_PATH}
setenv PYTHONPATH ${LAPACKDIR}/lib:${PYTHONPATH}
