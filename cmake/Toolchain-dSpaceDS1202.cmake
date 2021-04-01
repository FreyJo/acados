#
# (C) Copyright 2009 Johns Hopkins University (JHU), All Rights
# Reserved.
#
# --- begin cisst license - do not edit ---
# 
# This software is provided "as is" under an open source license, with
# no warranty.  The complete license can be found in license.txt and
# http://www.cisst.org/cisst/license.txt.
# 
# --- end cisst license ---

SET(CMAKE_SYSTEM_NAME QNX)
# set(CMAKE_SYSTEM_NAME dSpaceDS1202)
# list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

SET(CMAKE_SYSTEM_VERSION 6.5.0)
SET(CMAKE_SYSTEM_PROCESSOR ppc)
SET(TOOLCHAIN QNX)
SET(CMAKE_BUILD_TYPE STATIC_LIBRARY)

SET(CMAKE_SHARED_LIBRARY_PREFIX "lib")
SET(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
SET(CMAKE_STATIC_LIBRARY_PREFIX "lib")
SET(CMAKE_STATIC_LIBRARY_SUFFIX ".a")

IF(CMAKE_HOST_WIN32)
  SET(HOST_EXECUTABLE_SUFFIX ".exe")
ENDIF(CMAKE_HOST_WIN32)

FIND_PATH(QNX_HOST
  NAME usr/bin/make${HOST_EXECUTABLE_SUFFIX}
  #PATHS $ENV{QNX_HOST} C:/QNX650/host/win32/
  PATHS $ENV{QNX_HOST} C:/ProgramData/dSPACE/1184D92C-D928-4591-A1E9-B54339797C20/Compiler/QNX650_520/host/win32/x86/ #evtl hier Ende weg
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
)
  
FIND_PATH(QNX_TARGET
  NAME usr/include/qnx_errno.h
  #PATHS $ENV{QNX_TARGET} C:/QNX650/target/qnx6/
  PATHS $ENV{QNX_TARGET} C:/ProgramData/dSPACE/1184D92C-D928-4591-A1E9-B54339797C20/Compiler/QNX650_520/target/qnx6/
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
)

SET(ENV{QNX_HOST} ${QNX_HOST})
SET(ENV{QNX_TARGET} ${QNX_TARGET})

SET(CMAKE_MAKE_PROGRAM "${QNX_HOST}/usr/bin/make${HOST_EXECUTABLE_SUFFIX}"    CACHE PATH "QNX Make Program")
SET(CMAKE_SH           "${QNX_HOST}/usr/bin/sh${HOST_EXECUTABLE_SUFFIX}"      CACHE PATH "QNX shell Program")
SET(CMAKE_AR           "${QNX_HOST}/usr/bin/nto${CMAKE_SYSTEM_PROCESSOR}-ar${HOST_EXECUTABLE_SUFFIX}"      CACHE PATH "QNX ar Program")
SET(CMAKE_RANLIB       "${QNX_HOST}/usr/bin/nto${CMAKE_SYSTEM_PROCESSOR}-ranlib${HOST_EXECUTABLE_SUFFIX}"      CACHE PATH "QNX ranlib Program")
SET(CMAKE_NM           "${QNX_HOST}/usr/bin/nto${CMAKE_SYSTEM_PROCESSOR}-nm${HOST_EXECUTABLE_SUFFIX}"      CACHE PATH "QNX nm Program")
SET(CMAKE_OBJCOPY      "${QNX_HOST}/usr/bin/nto${CMAKE_SYSTEM_PROCESSOR}-objcopy${HOST_EXECUTABLE_SUFFIX}" CACHE PATH "QNX objcopy Program")
SET(CMAKE_OBJDUMP      "${QNX_HOST}/usr/bin/nto${CMAKE_SYSTEM_PROCESSOR}-objdump${HOST_EXECUTABLE_SUFFIX}" CACHE PATH "QNX objdump Program")
SET(CMAKE_LINKER       "${QNX_HOST}/usr/bin/nto${CMAKE_SYSTEM_PROCESSOR}-ld"     CACHE PATH "QNX Linker Program") #warum hier keinene hoste executalbe suffix? 
SET(CMAKE_STRIP        "${QNX_HOST}/usr/bin/nto${CMAKE_SYSTEM_PROCESSOR}-strip${HOST_EXECUTABLE_SUFFIX}"   CACHE PATH "QNX Strip Program")
#SET( CMAKE_LINKER       "${QNX_HOST}/usr/bin/qcc
#${HOST_EXECUTABLE_SUFFIX}"     CACHE PATH "QNX Linker Program" )

SET(CMAKE_C_COMPILER ${QNX_HOST}/usr/bin/nto${CMAKE_SYSTEM_PROCESSOR}-gcc${HOST_EXECUTABLE_SUFFIX})
SET(CMAKE_C_FLAGS_DEBUG "-g")
SET(CMAKE_C_FLAGS_MINSIZEREL "-Os -DNDEBUG")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

SET(CMAKE_CXX_COMPILER ${QNX_HOST}/usr/bin/nto${CMAKE_SYSTEM_PROCESSOR}-c++${HOST_EXECUTABLE_SUFFIX})
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

SET(CMAKE_FIND_ROOT_PATH ${QNX_TARGET}) 
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# acados flags
set(BLASFEO_TARGET "GENERIC" CACHE STRING "BLASFEO Target architecture")
set(HPIPM_TARGET "GENERIC" CACHE STRING "HPIPM Target architecture")
set(BUILD_SHARED_LIBS OFF CACHE STRING "Build shared libraries")
set(BLASFEO_EXAMPLES OFF CACHE BOOL "Examples disabled")
set(EXT_DEP OFF CACHE BOOL "Compile external dependencies in BLASFEO")
set(ACADOS_INSTALL_DIR "install" CACHE PATH  "Installation path to PROJECT_SOURCE_DIR")

# try to integrate dSpace Libraries
set(CMAKE_C_FLAGS "\"-I${DSPACE_RTLIB}\"")
set(CMAKE_INCLUDE_FLAG_C "-I")
set(CMAKE_INCLUDE_FLAG_CXX "-I")
# file(TO_CMAKE_PATH "C:\\Program Files\\dSPACE RCPHIL 2017-B" DSPACE_TOOLS)
set(DSPACE_RTLIB "C:/DualFuel/MLB/IdentificationOptimizationToolbox/acadosCrossCompile/acados/cmake/Platform/DS1202_RTLib")
# set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
# set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")

# add_definitions(-DWINDOWS_SKIP_PTR_ALIGNMENT_CHECK)
add_definitions(-DDSPACE_INCLUDES)
add_definitions(-D_DSHOST)
add_definitions(-D_DS1201)
add_definitions(-D_DS1202)
add_definitions(-DDS_PLATFORM_PPC)
add_definitions(-DDS_PLATFORM_SMARTRTK)
add_definitions(-DDS_PLATFORM_SMART)


# add_definitions(-D_INLINE)

