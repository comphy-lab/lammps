# Enable Fortran language support
enable_language(Fortran)

# Build arguments.
option(DOWNLOAD_RUNNER "Force download and build of RuNNer. If this is OFF\
  the arguments RUNNER_LIB_DIR, RUNNER_LIB_NAME and RUNNER_SHARED_LIB must\
  be set to meaningful values." ON
)
option(RUNNER_SHARED_LIB "Use pre-compiled shared/dynamic RuNNer library. Only \
  considered if DOWNLOAD_RUNNER is OFF." ON
)
set(RUNNER_LIB_DIR "" CACHE PATH "Directory containing \
  the RuNNer library. Only considered if DOWNLOAD_RUNNER is OFF."
)
set(RUNNER_LIB_NAME "libRuNNer_mpi" CACHE STRING "Name of the RuNNer library \
  (excluding file extension, this is controlled by RUNNER_SHARED_LIB). Only \
  considered if DOWNLOAD_RUNNER is OFF."
)

# FFT Library Selection
# Priority: 1. User choice via -D FFT=..., 2. MKL, 3. FFTW3, 4. KISSFFT.

# This block only runs if the user has NOT specified -D FFT=... on the command line.
if(NOT DEFINED FFT)
  # 1. Try to find MKL first.
  find_package(MKL QUIET)
  if(MKL_FOUND)
    set(FFT "MKL")
    message(STATUS "Auto-detected MKL. Using MKL for FFT by default.")
  else()
    # 2. If no MKL, try to find a standalone FFTW3.
    find_package(FFTW3 QUIET)
    if(FFTW3_FOUND)
      set(FFT "FFTW3")
      message(STATUS "Auto-detected FFTW3. Using standalone FFTW3 for FFT by default.")
    else()
      # 3. If nothing else is found, fall back to KISSFFT.
      set(FFT "KISS")
      message(STATUS "No MKL or FFTW3 found. Defaulting to KISS FFT.")
    endif()
  endif()
endif()

# Create the user-configurable cache variable with the determined default.
# This makes the option visible in tools like ccmake/cmake-gui and allows changes.
set(FFT_VALUES KISS FFTW3 MKL NVPL)
set(FFT ${FFT} CACHE STRING "FFT library for RUNNER package")
set_property(CACHE FFT PROPERTY STRINGS ${FFT_VALUES})

# Ensure the selected option is valid.
# validate_option(FFT FFT_VALUES) # Make sure this custom function is defined in your project
string(TOUPPER ${FFT} FFT_UPPER)
message(STATUS "Using ${FFT_UPPER} for FFT calculations.")

# Configure Project Based on Selected FFT Library

if(FFT_UPPER STREQUAL "MKL")
  find_package(MKL REQUIRED)
  target_compile_definitions(lammps PRIVATE -DFFT_MKL)
  option(FFT_MKL_THREADS "Use threaded MKL FFT" ON)
  if(FFT_MKL_THREADS)
    target_compile_definitions(lammps PRIVATE -DFFT_MKL_THREADS)
  endif()
  target_link_libraries(lammps PRIVATE MKL::MKL)

elseif(FFT_UPPER STREQUAL "FFTW3")
  find_package(FFTW3 REQUIRED)
  target_compile_definitions(lammps PRIVATE -DFFT_FFTW3)
  target_link_libraries(lammps PRIVATE FFTW3::FFTW3)

  # Check for OpenMP support in the found FFTW3 library
  if(FFTW3_OMP_LIBRARIES OR FFTW3F_OMP_LIBRARIES)
    option(FFT_FFTW_THREADS "Use threaded FFTW library" ON)
  else()
    option(FFT_FFTW_THREADS "Use threaded FFTW library" OFF)
  endif()

  if(FFT_FFTW_THREADS)
    if(FFTW3_OMP_LIBRARIES OR FFTW3F_OMP_LIBRARIES)
      target_compile_definitions(lammps PRIVATE -DFFT_FFTW_THREADS)
      target_link_libraries(lammps PRIVATE FFTW3::FFTW3_OMP)
    else()
      message(FATAL_ERROR "FFT_FFTW_THREADS is ON, but an OpenMP-enabled FFTW3 library was not found.")
    endif()
  endif()

elseif(FFT_UPPER STREQUAL "NVPL")
  find_package(nvpl_fft REQUIRED)
  target_compile_definitions(lammps PRIVATE -DFFT_NVPL)
  target_link_libraries(lammps PRIVATE nvpl::fftw)

else() # Fallback to KISSFFT
  if(NOT FFT_UPPER STREQUAL "KISS")
     message(WARNING "FFT option '${FFT}' not recognized. Falling back to KISSFFT.")
  endif()
  target_compile_definitions(lammps PRIVATE -DFFT_KISS)
endif()

if(DOWNLOAD_RUNNER)
  message(STATUS "DOWNLOAD_RUNNER is ON. Building RuNNer from source as a static library.")

  # Include ExternalProject module
  include(ExternalProject)

  # Force using the static library. RuNNer's cmake build always produces libRuNNer_mpi.a.
  if(BUILD_MPI)
    set(RUNNER_LIB_FULL_NAME "libRuNNer_mpi.a")
  else()
    set(RUNNER_LIB_FULL_NAME "libRuNNer.a")
  endif()

  # Add any custom CMake variables required by the RuNNer build system here.
  set(RUNNER_CMAKE_ARGS
    -DUSE_MPI=${BUILD_MPI}
    -DCMAKE_Fortran_FLAGS="-fPIC"
    -DENABLE_TESTS=no
  )

  ExternalProject_Add(runner_build
    GIT_REPOSITORY "git@gitlab.com:runner-suite/runner2.git"
    GIT_TAG "main"
    GIT_SHALLOW YES
    GIT_PROGRESS YES

    # Pass CMake arguments to RuNNer's build system
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/runner_install
      ${RUNNER_CMAKE_ARGS}

    # Define the build and install steps
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR>
    INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR>

    # Specify the location of the built library
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/runner_install/RuNNer/lib/${RUNNER_LIB_FULL_NAME}
  )

  # Create an IMPORTED library target for RuNNer
  add_library(RuNNer::RuNNer STATIC IMPORTED)
  set_target_properties(RuNNer::RuNNer PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/runner_install/RuNNer/lib/${RUNNER_LIB_FULL_NAME}"
  )

  # Add a dependency to ensure RuNNer is built before the main target
  add_dependencies(lammps runner_build)

else()
  message(STATUS "DOWNLOAD_RUNNER is OFF. Looking for a pre-compiled RuNNer.")

  # Check if the directory specified by RUNNER_LIB_DIR exists
  if(NOT IS_DIRECTORY ${RUNNER_LIB_DIR})
    message(FATAL_ERROR "The directory specified by RUNNER_LIB_DIR does not exist: ${RUNNER_LIB_DIR}")
  endif()

  # Determine the correct file extension based on whether the library is shared or not.
  if(RUNNER_SHARED_LIB)
    set(RUNNER_LIB_EXT ".so")
    set(RUNNER_LIB_TYPE SHARED)
  else()
    set(RUNNER_LIB_EXT ".a")
    set(RUNNER_LIB_TYPE STATIC)
  endif()

  set(RUNNER_LIB_FULL_NAME ${RUNNER_LIB_NAME}${RUNNER_LIB_EXT})
  message(STATUS "Looking for the ${RUNNER_LIB_TYPE} library: ${RUNNER_LIB_FULL_NAME}")

  # Find the RuNNer library in the specified path
  find_library(RuNNer_LIBRARY
    NAMES ${RUNNER_LIB_FULL_NAME}
    HINTS
      "${RUNNER_LIB_DIR}"       # 1st priority: User's manual cache variable
      ENV LD_LIBRARY_PATH       # 2nd priority: Search the shell's LD_LIBRARY_PATH
      ENV LIBRARY_PATH          # 3rd priority: Search the compiler's link path
    PATH_SUFFIXES lib lib64     # Automatically look in <hint>/lib if not found in <hint>
  )

  if(RuNNer_LIBRARY)
    message(STATUS "Found RuNNer library: ${RuNNer_LIBRARY}")

    # Create an IMPORTED library target for the found RuNNer library
    add_library(RuNNer::RuNNer ${RUNNER_LIB_TYPE} IMPORTED)
    set_target_properties(RuNNer::RuNNer PROPERTIES
      IMPORTED_LOCATION "${RuNNer_LIBRARY}"
    )
  else()
    message(FATAL_ERROR "Could not find the RuNNer library in the specified \
      RUNNER_LIB_DIR. Searched for ${RUNNER_LIB_FULL_NAME} at ${RUNNER_LIB_DIR}"
    )
  endif()

endif()

# Link lammps to the found RuNNer library
target_link_libraries(lammps PRIVATE RuNNer::RuNNer ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})