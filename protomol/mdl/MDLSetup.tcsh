#! /bin/tcsh

# USER: CHANGE THESE ACCORDINGLY
setenv PYTHONDIR /afs/nd.edu/user25/tcickovs/Research/PYTHON2.5
setenv LAPACKDIR /afs/nd.edu/user25/tcickovs/Research/LAPACK
setenv NUMPYDIR /afs/nd.edu/user25/tcickovs/Research/PYTHON2.5
setenv GNUPLOTDIR /afs/nd.edu/user25/tcickovs/Research/PYTHON2.5
setenv PYVERSION python2.5


setenv TMPPATH ${PWD}:${PWD}/interface/:${PWD}/interface/integrator:${PWD}/interface/integrator/base/:${PWD}/interface/integrator/hessian/:${PWD}/interface/integrator/normal:${PWD}/interface/integrator/leapfrog:${PWD}/interface/topology:${PWD}/interface/io/:${PWD}/interface/output:${PWD}/interface/force:${PWD}/interface/force/system:${PWD}/interface/force/bonded:${PWD}/interface/force/nonbonded:${PWD}/interface/type:${PWD}/interface/base:${PWD}/src/:${PWD}/src/factories:${PWD}/src/modifiers:${PWD}/src/propagators:${PWD}/src/propagators/objects:${PWD}/src/toplevel:${GNUPLOTDIR}/lib/${PYVERSION}/site-packages/Gnuplot:${NUMPYDIR}/lib/${PYVERSION}/site-packages/:${NUMPYDIR}/lib/${PYVERSION}/site-packages/numpy/

if ($?PYTHONPATH != 0) then
   setenv PYTHONPATH ${TMPPATH}:${PYTHONPATH}
else
   setenv PYTHONPATH ${TMPPATH}
endif

setenv MDLROOT ${PWD}
setenv PATH ${PYTHONDIR}/bin:${PATH}
setenv LD_LIBRARY_PATH ${PYTHONDIR}/lib/${PYVERSION}:${LD_LIBRARY_PATH}
setenv PYTHONPATH ${PYTHONDIR}/lib/${PYVERSION}:${PYTHONPATH}
setenv LD_LIBRARY_PATH ${PWD}/interface/:${LAPACKDIR}/lib:${LD_LIBRARY_PATH}
setenv PYTHONPATH ${LAPACKDIR}/lib:${PYTHONPATH}
