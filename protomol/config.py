if env['PLATFORM'] == 'win32':
    env.Append(LIBS = 'wsock32')
else:
    env.Append(LIBS = 'pthread')

if env['CC'] == 'gcc':
    if os.environ.has_key('ATLAS_HOME'):
        env.Append(LIBPATH = [os.environ['ATLAS_HOME']])
    else:
        env.Append(LIBPATH = ['/usr/lib/atlas'])

if globals().has_key('use_openmm') and use_openmm:
    if os.environ.has_key('OPENMM_HOME') == False:
        print 'OPENMM_HOME is not set'
    else:
        have_openmm = 0

        openmm_home = os.environ['OPENMM_HOME']
        env.Append(CPPPATH = [openmm_home + os.sep +'include'])
        env.Append(LIBPATH = [openmm_home + os.sep + 'lib'   ])

        if (globals().has_key('use_openmm_reference') and use_openmm_reference) or (globals().has_key('use_openmm_cuda') and use_openmm_cuda):
            if conf.CheckLib('OpenMM_d'):
                env.Append(CPPDEFINES = ['HAVE_OPENMM'])
                have_openmm = 1

        if globals().has_key('use_openmm_cuda') and use_openmm_cuda:
            if os.environ.has_key('CUDA_HOME') == False:
                print 'CUDA_HOME is not set'
            else:
                openmm_home = os.environ['CUDA_HOME']
                env.Append(LIBPATH = [openmm_home + os.sep + 'lib'])
            
                if conf.CheckLib('cudart') and conf.CheckLib('OpenMMCuda_d'):
                    env.Append(CPPDEFINES = ['HAVE_OPENMM'])
                    have_openmm = 1

# LAPACK
if globals().has_key('use_lapack') and use_lapack:
    if os.environ.has_key('LAPACK_HOME'):
        env.Append(CPPPATH = [os.environ['LAPACK_HOME']])
        env.Append(LIBPATH = [os.environ['LAPACK_HOME']])
  
    have_lapack = 0
    if conf.CheckLib('lapack'):
        env.Append(CPPDEFINES = ['HAVE_LAPACK'])
        have_lapack = 1


    if env['CC'] == 'gcc' and have_lapack:
        # BLAS
        if os.environ.has_key('BLAS_HOME'):
            env.Append(LIBPATH = [os.environ['BLAS_HOME']])

        have_blas = 0  
        if conf.CheckLib('blas'):
            have_blas = 1


        # G2C
        if os.environ.has_key('G2C_HOME'):
            env.Append(LIBPATH = [os.environ['G2C_HOME']])

        have_g2c = 0  
        if conf.CheckLib('g2c'):
            have_g2c = 1


if globals().has_key('use_simtk_lapack') and use_simtk_lapack:
    # SimTK LAPACK
    if os.environ.has_key('SIMTK_LAPACK_HOME'):
        env.Append(CPPPATH = [os.environ['SIMTK_LAPACK_HOME'] + '/include'])
        env.Append(LIBPATH = [os.environ['SIMTK_LAPACK_HOME'] + '/lib'])
  
    have_simtk_lapack = 0
    if conf.CheckLib('SimTKlapack') and conf.CheckCXXHeader('SimTKlapack.h'):
        env.Append(CPPDEFINES = ['HAVE_SIMTK_LAPACK'])
        have_simtk_lapack = 1
