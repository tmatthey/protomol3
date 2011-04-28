cmake_minimum_required( VERSION 2.6 )

project( ProtoMol )

set( CMAKE_VERBOSE_MAKEFILE Off )

option( BUILD_FAH "Build for Folding@Home" Off )
if( BUILD_FAH )	
endif()

option( BUILD_LAPACK "Build with Lapack" On )
if( BUILD_LAPACK )	
	set( LAPACK_LIBS "" )
	set( LAPACK_DEFS "" )
	
	set( BUILD_LAPACK_LIBDIR "" CACHE PATH "Lapack type library directory" )
	set( BUILD_LAPACK_INCDIR "" CACHE PATH "Lapack type include directory" )
	set( BUILD_LAPACK_TYPE "lapack" CACHE STRING "Lapack type to compile against. Choices are: lapack simtk mkl" )
		
	if( "${BUILD_LAPACK_INCDIR}" STRGREATER "" )
		include_directories( "${BUILD_LAPACK_INCDIR}" )
	endif()
	
	if( ${BUILD_LAPACK_TYPE} MATCHES "simtk" )
		
		find_library( ZLIB_LAPACK "SimTKlapack" HINTS "${BUILD_LAPACK_LIBDIR}" )
		
		set( LAPACK_LIBS ${ZLIB_LAPACK} )
		add_definitions( "/DHAVE_SIMTK_LAPACK" )
		#set( LAPACK_DEFS "/DHAVE_SIMTK_LAPACK" )
	else()
		if( ${BUILD_LAPACK_TYPE} MATCHES "mkl" )		
			find_library( ZLIB_MKLCORE "mkl_core" HINTS "${BUILD_BLAS_LIBDIR}" )
			find_library( ZLIB_MKLSEQ "mkl_sequential" HINTS "${BUILD_LAPACK_LIBDIR}" )
			find_library( ZLIB_MKLINTEL "mkl_intel_c" HINTS "${BUILD_LAPACK_LIBDIR}" )
			
			list( APPEND LAPACK_LIBS ${ZLIB_MKLINTEL} )
			list( APPEND LAPACK_LIBS ${ZLIB_MKLSEQ} )
			list( APPEND LAPACK_LIBS ${ZLIB_MKLCORE} )
			add_definitions( "/DHAVE_MKL_LAPACK" )
			#set( LAPACK_DEFS "/DHAVE_MKL_LAPACK" )
		else()
			set( BUILD_BLAS_LIBDIR "" CACHE PATH "Blas library directory" )
			set( BUILD_BLAS_INCDIR "" CACHE PATH "Blas include directory" )
			
			if( "${BUILD_BLAS_INCDIR}" STRGREATER "" )
				include_directories( "${BUILD_BLAS_INCDIR}" )
			endif()
		
			find_library( ZLIB_BLAS "blas" HINTS "${BUILD_BLAS_LIBDIR}" )
			find_library( ZLIB_LAPACK "lapack" HINTS "${BUILD_LAPACK_LIBDIR}" )
			
			list( APPEND LAPACK_LIBS ${ZLIB_BLAS} )
			list( APPEND LAPACK_LIBS ${ZLIB_LAPACK} )
			add_definitions( "/DHAVE_LAPACK" )
			#set( LAPACK_DEFS "/DHAVE_LAPACK" )
		endif() 
	endif()
endif()

# Compilation
include_directories( src )

file( GLOB_RECURSE LIBRARY_SOURCES "src/*.cpp")
file( GLOB_RECURSE LIBRARY_HEADERS "src/*.h" )

file( GLOB_RECURSE DEBUG_SOURCES "src/protomol/debug/*.cpp" )

foreach( f ${DEBUG_SOURCES} )
	list( REMOVE_ITEM LIBRARY_SOURCES ${f} )
endforeach()

list( REMOVE_ITEM LIBRARY_SOURCES "src/protomol/main.cpp" )
list( REMOVE_ITEM LIBRARY_SOURCES "src/protomol/modules.cpp" )

add_definitions( ${LAPACK_DEFS} )
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LAPACK_DEFS} /O3 /QaxAVX /arch:IA32" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LAPACK_DEFS} /O3 /QaxAVX /arch:IA32" )

add_library( protomol STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS} )
target_link_libraries( protomol ${LAPACK_LIBS} )
set_target_properties( protomol PROPERTIES COMPILE_FLAGS "${LAPACK_DEFS} /O3 /QaxAVX /arch:IA32" )

add_executable( ProtoMolMain "src/protomol/main.cpp" "src/protomol/modules.cpp" )
target_link_libraries( ProtoMolMain protomol )
set_target_properties( ProtoMolMain PROPERTIES COMPILE_FLAGS "${LAPACK_DEFS} /O3 /QaxAVX /arch:IA32" )