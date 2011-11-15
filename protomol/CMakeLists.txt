project( ProtoMol )
cmake_minimum_required( VERSION 2.6 )

option( BUILD_FAH "Build for Folding@Home" Off )
if( BUILD_FAH )
endif()

set( LIBS "" )

option( BUILD_LAPACK "Build with Lapack" On )
if( BUILD_LAPACK )
	set( BUILD_LAPACK_LIBDIR "" CACHE PATH "Lapack type library directory" )
	set( BUILD_LAPACK_INCDIR "" CACHE PATH "Lapack type include directory" )
	set( BUILD_LAPACK_TYPE "lapack" CACHE STRING "Lapack type to compile against. Choices are: lapack simtk mkl" )

	if( "${BUILD_LAPACK_INCDIR}" STRGREATER "" )
		include_directories( "${BUILD_LAPACK_INCDIR}" )
	endif()

	if( ${BUILD_LAPACK_TYPE} MATCHES "simtk" )
		find_library( ZLIB_LAPACK "SimTKlapack" HINTS "${BUILD_LAPACK_LIBDIR}" )

		set( LIBS ${LIBS} ${ZLIB_LAPACK} )
		add_definitions( "-DHAVE_SIMTK_LAPACK" )
	else()
		if( ${BUILD_LAPACK_TYPE} MATCHES "mkl" )
			find_library( ZLIB_MKLCORE "mkl_core" HINTS "${BUILD_BLAS_LIBDIR}" )
			find_library( ZLIB_MKLSEQ "mkl_sequential" HINTS "${BUILD_LAPACK_LIBDIR}" )
			find_library( ZLIB_MKLINTEL "mkl_intel_c" HINTS "${BUILD_LAPACK_LIBDIR}" )

			set( LIBS ${LIBS} ${ZLIB_MKLINTEL} ${ZLIB_MKLSEQ} ${ZLIB_MKLCORE} )
			add_definitions( "-DHAVE_MKL_LAPACK" )
		else()
			set( BUILD_BLAS_LIBDIR "" CACHE PATH "Blas library directory" )
			set( BUILD_BLAS_INCDIR "" CACHE PATH "Blas include directory" )

			if( "${BUILD_BLAS_INCDIR}" STRGREATER "" )
				include_directories( "${BUILD_BLAS_INCDIR}" )
			endif()

			find_library( ZLIB_BLAS "blas" HINTS "${BUILD_BLAS_LIBDIR}" )
			find_library( ZLIB_LAPACK "lapack" HINTS "${BUILD_LAPACK_LIBDIR}" )

			set( LIBS ${LIBS} ${ZLIB_BLAS} ${ZLIB_LAPACK} )
			add_definitions( "-DHAVE_LAPACK" )
		endif()
	endif()
endif()

option( BUILD_OPENMM "Build with OpenMM" On )
if( BUILD_OPENMM )
	add_definitions( "-DHAVE_OPENMM" )
	
	set( OPENMM_DIR "" CACHE PATH "OpenMM Install Directory")
	if( OPENMM_DIR )
		link_directories( "${OPENMM_DIR}/lib/" )
		include_directories( "${OPENMM_DIR}/include/" )
	endif()

	find_library( OPENMM_LIB "OpenMM" HINT ${OPENMM_DIR}/lib )
	set( LIBS ${LIBS} ${OPENMM_LIB} )
	
	set( NMLPLUGIN_DIR "" CACHE PATH "NML OpenMM Install Directory")
	if( NMLPLUGIN_DIR )
		include_directories( "${NMLPLUGIN_DIR}/include/" )
	endif()
endif()

# Compilation
include_directories( src )

# Library Headers
file( GLOB HEADERS "src/protomol/*.h" )
source_group( "include" FILES ${HEADERS} )
file( GLOB BASE_HEADERS "src/protomol/base/*.h" )
source_group( "include\\Base" FILES ${BASE_HEADERS} )
file( GLOB CONFIG_HEADERS "src/protomol/config/*.h" )
source_group( "include\\Config" FILES ${CONFIG_HEADERS} )
file( GLOB FACTORY_HEADERS "src/protomol/factory/*.h" )
source_group( "include\\Factory" FILES ${FACTORY_HEADERS} )
file( GLOB MODIFIER_HEADERS "src/protomol/modifier/*.h" )
source_group( "include\\Modifier" FILES ${MODIFIER_HEADERS} )
file( GLOB MODULE_HEADERS "src/protomol/module/*.h" )
source_group( "include\\Module" FILES ${MODULE_HEADERS} )
file( GLOB OS_HEADERS "src/protomol/os/*.h" )
source_group( "include\\OS" FILES ${OS_HEADERS} )
file( GLOB OUTPUT_HEADERS "src/protomol/output/*.h" )
source_group( "include\\Output" FILES ${OUTPUT_HEADERS} )
file( GLOB PARALLEL_HEADERS "src/protomol/parallel/*.h" )
source_group( "include\\Parallel" FILES ${PARALLEL_HEADERS} )
file( GLOB SWITCH_HEADERS "src/protomol/switch/*.h" )
source_group( "include\\Switch" FILES ${SWITCH_HEADERS} )
file( GLOB TOPOLOGY_HEADERS "src/protomol/topology/*.h" )
source_group( "include\\Topology" FILES ${TOPOLOGY_HEADERS} )
file( GLOB TYPE_HEADERS "src/protomol/type/*.h" )
source_group( "include\\Type" FILES ${TYPE_HEADERS} )

set( LIBRARY_BASE_HEADERS ${HEADERS} ${BASE_HEADERS} ${CONFIG_HEADERS}
	${FACTORY_HEADERS} ${MODIFIER_HEADERS} ${MODULE_HEADERS} ${OS_HEADERS}
	${OUTPUT_HEADERS} ${PARALLEL_HEADERS} ${SWITCH_HEADERS} ${TOPOLOGY_HEADERS}
	${TYPE_HEADERS} )
	
if( NOT MSVC )
	file( GLOB DEBUG_HEADERS "src/protomol/debug/*.h" )
	source_group( "include\\Debug" FILES ${DEBUG_HEADERS} )

	set( LIBRARY_BASE_HEADERS ${LIBRARY_BASE_HEADERS} ${DEBUG_HEADERS} )
endif( NOT MSVC )

# Library Force Headers
file( GLOB FORCE_HEADERS "src/protomol/force/*.h" )
source_group( "include\\Force" FILES ${FORCE_HEADERS} )
file( GLOB BONDED_FORCE_HEADERS "src/protomol/force/bonded/*.h" )
source_group( "include\\Force\\Bonded" FILES ${BONDED_FORCE_HEADERS} )
file( GLOB BORN_FORCE_HEADERS "src/protomol/force/born/*.h" )
source_group( "include\\Force\\Born" FILES ${BORN_FORCE_HEADERS} )
file( GLOB COULOMB_FORCE_HEADERS "src/protomol/force/coulomb/*.h" )
source_group( "include\\Force\\Coulomb" FILES ${COULOMB_FORCE_HEADERS} )
file( GLOB EXTENDED_FORCE_HEADERS "src/protomol/force/extended/*.h" )
source_group( "include\\Force\\Extended" FILES ${EXTENDED_FORCE_HEADERS} )
file( GLOB GB_FORCE_HEADERS "src/protomol/force/gb/*.h" )
source_group( "include\\Force\\GB" FILES ${GB_FORCE_HEADERS} )
file( GLOB HESSIAN_FORCE_HEADERS "src/protomol/force/hessian/*.h" )
source_group( "include\\Force\\Hessian" FILES ${HESSIAN_FORCE_HEADERS} )
file( GLOB NONBONDED_FORCE_HEADERS "src/protomol/force/nonbonded/*.h" )
source_group( "include\\Force\\NonBonded" FILES ${NONBONDED_FORCE_HEADERS} )
file( GLOB SYSTEM_FORCE_HEADERS "src/protomol/force/system/*.h" )
source_group( "include\\Force\\System" FILES ${SYSTEM_FORCE_HEADERS} )
file( GLOB TABLE_FORCE_HEADERS "src/protomol/force/table/*.h" )
source_group( "include\\Force\\Table" FILES ${TABLE_FORCE_HEADERS} )

set( LIBRARY_FOURCE_HEADERS ${FORCE_HEADERS} ${BONDED_FORCE_HEADERS}
	${BORN_FORCE_HEADERS} ${COULOMB_FORCE_HEADERS} ${EXTENDED_FORCE_HEADERS}
	${GB_FORCE_HEADERS} ${HESSIAN_FORCE_HEADERS} ${NONBONDED_FORCE_HEADERS}
	${SYSTEM_FORCE_HEADERS} ${TABLE_FORCE_HEADERS} )

# Library Integrator Headers
file( GLOB INTEGRATOR_HEADERS "src/protomol/integrator/*.h" )
source_group( "include\\Integrator" FILES ${INTEGRATOR_HEADERS} )
file( GLOB BASE_INTEGRATOR_HEADERS "src/protomol/integrator/base/*.h" )
source_group( "include\\Integrator\\Base" FILES ${BASE_INTEGRATOR_HEADERS} )
file( GLOB HESSIAN_INTEGRATOR_HEADERS "src/protomol/integrator/hessian/*.h" )
source_group( "include\\Integrator\\Hessian" FILES ${HESSIAN_INTEGRATOR_HEADERS} )
file( GLOB LEAPFROG_INTEGRATOR_HEADERS "src/protomol/integrator/leapfrog/*.h" )
source_group( "include\\Integrator\\Leapfrog" FILES ${LEAPFROG_INTEGRATOR_HEADERS} )
file( GLOB NORMAL_INTEGRATOR_HEADERS "src/protomol/integrator/normal/*.h" )
source_group( "include\\Integrator\\Normal" FILES ${NORMAL_INTEGRATOR_HEADERS} )
file( GLOB OPENMM_INTEGRATOR_HEADERS "src/protomol/integrator/openMM/*.h" )
source_group( "include\\Integrator\\OpenMM" FILES ${OPENMM_INTEGRATOR_HEADERS} )

set( LIBRARY_INTEGRATOR_HEADERS ${INTEGRATOR_HEADERS} ${BASE_INTEGRATOR_HEADERS}
	${HESSIAN_INTEGRATOR_HEADERS} ${LEAPFROG_INTEGRATOR_HEADERS} 
	${NORMAL_INTEGRATOR_HEADERS} )
	
if( BUILD_OPENMM )
	set( LIBRARY_INTEGRATOR_HEADERS ${LIBRARY_INTEGRATOR_HEADERS} ${OPENMM_INTEGRATOR_HEADERS} )
endif( BUILD_OPENMM )

# Library IO Headers
file( GLOB IO_HEADERS "src/protomol/io/*.h" )
source_group( "include\\IO" FILES ${IO_HEADERS} )
file( GLOB GROMACS_IO_HEADERS "src/protomol/io/gromacs/*.h" )
source_group( "include\\IO\\Gromacs" FILES ${GROMACS_IO_HEADERS} )

set( LIBRARY_IO_HEADERS ${IO_HEADERS} ${GROMACS_IO_HEADERS} )

# Combined Headers
set( LIBRARY_HEADERS ${LIBRARY_BASE_HEADERS} ${LIBRARY_FOURCE_HEADERS} ${LIBRARY_INTEGRATOR_HEADERS} ${LIBRARY_IO_HEADERS} )

# Library sources
file( GLOB SOURCES "src/protomol/*.cpp" )
source_group( "src" FILES ${SOURCES} )
file( GLOB BASE_SOURCES "src/protomol/base/*.cpp" )
source_group( "src\\Base" FILES ${BASE_SOURCES} )
file( GLOB CONFIG_SOURCES "src/protomol/config/*.cpp" )
source_group( "src\\Config" FILES ${CONFIG_SOURCES} )
file( GLOB FACTORY_SOURCES "src/protomol/factory/*.cpp" )
source_group( "src\\Factory" FILES ${FACTORY_SOURCES} )
file( GLOB MODIFIER_SOURCES "src/protomol/modifier/*.cpp" )
source_group( "src\\Modifier" FILES ${MODIFIER_SOURCES} )
file( GLOB MODULE_SOURCES "src/protomol/module/*.cpp" )
source_group( "src\\Module" FILES ${MODULE_SOURCES} )
file( GLOB OS_SOURCES "src/protomol/os/*.cpp" )
source_group( "src\\OS" FILES ${OS_SOURCES} )
file( GLOB OUTPUT_SOURCES "src/protomol/output/*.cpp" )
source_group( "src\\Output" FILES ${OUTPUT_SOURCES} )
file( GLOB PARALLEL_SOURCES "src/protomol/parallel/*.cpp" )
source_group( "src\\Parallel" FILES ${PARALLEL_SOURCES} )
file( GLOB SWITCH_SOURCES "src/protomol/switch/*.cpp" )
source_group( "src\\Switch" FILES ${SWITCH_SOURCES} )
file( GLOB TOPOLOGY_SOURCES "src/protomol/topology/*.cpp" )
source_group( "src\\Topology" FILES ${TOPOLOGY_SOURCES} )
file( GLOB TYPE_SOURCES "src/protomol/type/*.cpp" )
source_group( "src\\Type" FILES ${TYPE_SOURCES} )

set( LIBRARY_BASE_SOURCES ${SOURCES} ${BASE_SOURCES} ${CONFIG_SOURCES}
	${FACTORY_SOURCES} ${MODIFIER_SOURCES} ${MODULE_SOURCES} ${OS_SOURCES}
	${OUTPUT_SOURCES} ${PARALLEL_SOURCES} ${SWITCH_SOURCES} ${TOPOLOGY_SOURCES}
	${TYPE_SOURCES} )

if( NOT MSVC )
	file( GLOB DEBUG_SOURCES "src/protomol/debug/*.cpp" )
	source_group( "src\\Debug" FILES ${DEBUG_SOURCES} )

	set( LIBRARY_BASE_SOURCES ${LIBRARY_BASE_SOURCES} ${DEBUG_SOURCES} )
endif( NOT MSVC )

# Library Force sources
file( GLOB FORCE_SOURCES "src/protomol/force/*.cpp" )
source_group( "src\\Force" FILES ${FORCE_SOURCES} )
file( GLOB BONDED_FORCE_SOURCES "src/protomol/force/bonded/*.cpp" )
source_group( "src\\Force\\Bonded" FILES ${BONDED_FORCE_SOURCES} )
file( GLOB BORN_FORCE_SOURCES "src/protomol/force/born/*.cpp" )
source_group( "src\\Force\\Born" FILES ${BORN_FORCE_SOURCES} )
file( GLOB COULOMB_FORCE_SOURCES "src/protomol/force/coulomb/*.cpp" )
source_group( "src\\Force\\Coulomb" FILES ${COULOMB_FORCE_SOURCES} )
file( GLOB EXTENDED_FORCE_SOURCES "src/protomol/force/extended/*.cpp" )
source_group( "src\\Force\\Extended" FILES ${EXTENDED_FORCE_SOURCES} )
file( GLOB GB_FORCE_SOURCES "src/protomol/force/gb/*.cpp" )
source_group( "src\\Force\\GB" FILES ${GB_FORCE_SOURCES} )
file( GLOB HESSIAN_FORCE_SOURCES "src/protomol/force/hessian/*.cpp" )
source_group( "src\\Force\\Hessian" FILES ${HESSIAN_FORCE_SOURCES} )
file( GLOB NONBONDED_FORCE_SOURCES "src/protomol/force/nonbonded/*.cpp" )
source_group( "src\\Force\\NonBonded" FILES ${NONBONDED_FORCE_SOURCES} )
file( GLOB SYSTEM_FORCE_SOURCES "src/protomol/force/system/*.cpp" )
source_group( "src\\Force\\System" FILES ${SYSTEM_FORCE_SOURCES} )
file( GLOB TABLE_FORCE_SOURCES "src/protomol/force/table/*.cpp" )
source_group( "src\\Force\\Table" FILES ${TABLE_FORCE_SOURCES} )

set( LIBRARY_FOURCE_SOURCES ${FORCE_SOURCES} ${BONDED_FORCE_SOURCES}
	${BORN_FORCE_SOURCES} ${COULOMB_FORCE_SOURCES} ${EXTENDED_FORCE_SOURCES}
	${GB_FORCE_SOURCES} ${HESSIAN_FORCE_SOURCES} ${NONBONDED_FORCE_SOURCES}
	${SYSTEM_FORCE_SOURCES} ${TABLE_FORCE_SOURCES} )

# Library Integrator sources
file( GLOB INTEGRATOR_SOURCES "src/protomol/integrator/*.cpp" )
source_group( "src\\Integrator" FILES ${INTEGRATOR_SOURCES} )
file( GLOB BASE_INTEGRATOR_SOURCES "src/protomol/integrator/base/*.cpp" )
source_group( "src\\Integrator\\Base" FILES ${BASE_INTEGRATOR_SOURCES} )
file( GLOB HESSIAN_INTEGRATOR_SOURCES "src/protomol/integrator/hessian/*.cpp" )
source_group( "src\\Integrator\\Hessian" FILES ${HESSIAN_INTEGRATOR_SOURCES} )
file( GLOB LEAPFROG_INTEGRATOR_SOURCES "src/protomol/integrator/leapfrog/*.cpp" )
source_group( "src\\Integrator\\Leapfrog" FILES ${LEAPFROG_INTEGRATOR_SOURCES} )
file( GLOB NORMAL_INTEGRATOR_SOURCES "src/protomol/integrator/normal/*.cpp" )
source_group( "src\\Integrator\\Normal" FILES ${NORMAL_INTEGRATOR_SOURCES} )
file( GLOB OPENMM_INTEGRATOR_SOURCES "src/protomol/integrator/openMM/*.cpp" )
source_group( "src\\Integrator\\OpenMM" FILES ${OPENMM_INTEGRATOR_SOURCES} )

set( LIBRARY_INTEGRATOR_SOURCES ${INTEGRATOR_SOURCES} ${BASE_INTEGRATOR_SOURCES}
	${HESSIAN_INTEGRATOR_SOURCES} ${LEAPFROG_INTEGRATOR_SOURCES} 
	${NORMAL_INTEGRATOR_SOURCES} )

if( BUILD_OPENMM )
	set( LIBRARY_INTEGRATOR_SOURCES ${LIBRARY_INTEGRATOR_SOURCES} ${OPENMM_INTEGRATOR_SOURCES} )
endif( BUILD_OPENMM )

# Library IO sources
file( GLOB IO_SOURCES "src/protomol/io/*.cpp" )
source_group( "src\\IO" FILES ${IO_SOURCES} )
file( GLOB GROMACS_IO_SOURCES "src/protomol/io/gromacs/*.cpp" )
source_group( "src\\IO\\Gromacs" FILES ${GROMACS_IO_SOURCES} )

set( LIBRARY_IO_SOURCES ${IO_SOURCES} ${GROMACS_IO_SOURCES} )

# Combined sources
set( LIBRARY_SOURCES ${LIBRARY_BASE_SOURCES} ${LIBRARY_FOURCE_SOURCES} ${LIBRARY_INTEGRATOR_SOURCES} ${LIBRARY_IO_SOURCES} )

# Library Compilation
add_library( protomol STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS} )
target_link_libraries( protomol ${LIBS} )

# Application Compilation
add_executable( ProtoMol "src/protomol/main.cpp" "src/protomol/modules.cpp" )
target_link_libraries( ProtoMol protomol )
