# Try to locate CVMix's installation prefix.
find_path(CVMix_PREFIX
    NAMES src/cvmix_driver.F90
    HINTS "$ENV{CVMIX_ROOT}"
    DOC "Installation prefix for the Community Vertical Mixing Project - cvmix.github.io"
)

# Find CVMix library
find_library(CVMix_LIBRARIES
    NAMES cvmix
    HINTS ${CVMix_PREFIX}/lib
    DOC "CVMix libraries")

# Store configurable path of CVMix include directory
find_path(CVMix_INCLUDE_DIRS
    NAME cvmix_kpp.mod
    HINTS ${CVMix_PREFIX}/include
    DOC "CVMix include directory"
)

mark_as_advanced(CVMix_LIBRARIES CVMix_INCLUDE_DIRS)

# For backward compatibility:
set(CVMix_LIBRARY CVMix_LIBRARIES)
set(CVMix_INCLUDE_DIR CVMix_INCLUDE_DIRS)
