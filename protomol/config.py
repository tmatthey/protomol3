if env['CC'] == 'gcc':
    if os.environ.has_key('ATLAS_HOME'):
        env.Append(LIBPATH = [os.environ['ATLAS_HOME']])
    else:
        env.Append(LIBPATH = ['/usr/lib/atlas'])


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
