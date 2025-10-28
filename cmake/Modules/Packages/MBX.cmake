# MBX v1.3.0+ support for MBX package
# Based off of PLUMED.cmake

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

# for supporting multiple concurrent mbx installations for debugging and testing
set(MBX_SUFFIX "" CACHE STRING "Suffix for MBX library")
mark_as_advanced(MBX_SUFFIX)

set(MBX_CONFIG_CC  ${CMAKE_C_COMPILER})
set(MBX_CONFIG_CXX  ${CMAKE_CXX_COMPILER})
if(BUILD_MPI)
  set(MBX_CONFIG_MPI "--enable-mpi")
else()
  set(MBX_CONFIG_MPI "--disable-mpi")
endif()


# TODO: change to static release tarball
set(MBX_URL "https://github.com/paesanilab/MBX/archive/refs/tags/v1.3.0.tar.gz"
  CACHE STRING "URL for MBX tarball")
set(MBX_MD5 "7cfbf221f9c249c364c7f47769d0d768" CACHE STRING "MD5 checksum of MBX tarball")

mark_as_advanced(MBX_URL)
mark_as_advanced(MBX_MD5)



set(MBX_MODE "static" CACHE STRING "Linkage mode for MBX library")
set(MBX_MODE_VALUES static shared runtime)
set_property(CACHE MBX_MODE PROPERTY STRINGS ${MBX_MODE_VALUES})
validate_option(MBX_MODE MBX_MODE_VALUES)
string(TOUPPER ${MBX_MODE} MBX_MODE)

set(MBX_LINK_LIBS)
if(MBX_MODE STREQUAL "STATIC")
  find_package(FFTW3 QUIET)
  if(FFTW3_FOUND)
    list(APPEND MBX_LINK_LIBS FFTW3::FFTW3)
  endif()
endif()

find_package(PkgConfig QUIET)
set(DOWNLOAD_MBX_DEFAULT ON)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(MBX QUIET mbx)
  if(MBX_FOUND)
    set(DOWNLOAD_MBX_DEFAULT OFF)
  endif()
endif()

option(DOWNLOAD_MBX "Download MBX package instead of using an already installed one" ${DOWNLOAD_MBX_DEFAULT})
if(DOWNLOAD_MBX)
  message(STATUS "MBX download requested - we will build our own")
  set(MBX_BUILD_BYPRODUCTS "<INSTALL_DIR>/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mbx${CMAKE_STATIC_LIBRARY_SUFFIX}")

  # MBX perform autoreconf -fi

  message(STATUS "MBX_CONFIG_MPI: ${MBX_CONFIG_MPI}")
  message(MBX_CONFIG_CXX: ${MBX_CONFIG_CXX})
  message(CPPFLAGS: ${MBX_CONFIG_CPPFLAGS})

  include(ExternalProject)
  ExternalProject_Add(mbx_build
    URL     ${MBX_URL}
    URL_MD5 ${MBX_MD5}
    CONFIGURE_COMMAND --prefix=<INSTALL_DIR>
                      ${MBX_CONFIG_MPI}
                      CXX=${MBX_CONFIG_CXX}
                      CC=${MBX_CONFIG_CC}
                      # CPPFLAGS=${MBX_CONFIG_CPPFLAGS}
    BUILD_BYPRODUCTS ${MBX_BUILD_BYPRODUCTS}
  )
  ExternalProject_get_property(mbx_build INSTALL_DIR)
  add_library(LAMMPS::MBX UNKNOWN IMPORTED)
  add_dependencies(LAMMPS::MBX mbx_build)
  set_target_properties(LAMMPS::MBX PROPERTIES
    IMPORTED_LOCATION ${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mbx${CMAKE_STATIC_LIBRARY_SUFFIX}
    INTERFACE_LINK_LIBRARIES "${MBX_LINK_LIBS};${CMAKE_DL_LIBS}"
  )

  set_target_properties(LAMMPS::MBX PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
  if(CMAKE_PROJECT_NAME STREQUAL "lammps")
    target_link_libraries(lammps PRIVATE LAMMPS::MBX)
  endif()
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(MBX REQUIRED mbx${MBX_SUFFIX})
  add_library(LAMMPS::MBX INTERFACE IMPORTED)

  file(MAKE_DIRECTORY ${MBX_INCLUDE_DIRS})

  target_include_directories(LAMMPS::MBX INTERFACE ${MBX_INCLUDE_DIRS})
  target_link_directories(LAMMPS::MBX INTERFACE ${MBX_LIBRARY_DIRS})
  target_link_libraries(LAMMPS::MBX INTERFACE "${MBX_LINK_LIBS};${MBX_LIBRARIES};${CMAKE_DL_LIBS}")

  if(CMAKE_PROJECT_NAME STREQUAL "lammps")
    target_link_libraries(lammps PUBLIC LAMMPS::MBX)
  endif()
endif()

