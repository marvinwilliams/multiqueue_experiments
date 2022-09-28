#[=======================================================================[.rst:
FindPAPI
-------

Finds the PAPI library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``PAPI::PAPI``

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``PAPI_FOUND``
  True if the system has the PAPI library.
``PAPI_INCLUDE_DIRS``
  Include directories needed to use PAPI.
``PAPI_LIBRARIES``
  Libraries needed to link to PAPI.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``PAPI_INCLUDE_DIR``
  The directory containing ``papi.h``.
``PAPI_LIBRARY``
  The path to the PAPI library.

#]=======================================================================]

find_path(PAPI_INCLUDE_DIR papi.h)
find_library(PAPI_LIBRARY papi)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  PAPI
  FOUND_VAR PAPI_FOUND
  REQUIRED_VARS
  PAPI_LIBRARY
  PAPI_INCLUDE_DIR
)

if(PAPI_FOUND)
  set(PAPI_LIBRARIES ${PAPI_LIBRARY})
  set(PAPI_INCLUDE_DIRS ${PAPI_INCLUDE_DIR})
  if(NOT TARGET PAPI::PAPI)
    add_library(PAPI::PAPI UNKNOWN IMPORTED)
    set_target_properties(
      PAPI::PAPI PROPERTIES
      IMPORTED_LOCATION "${PAPI_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${PAPI_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(
  PAPI_INCLUDE_DIR
  PAPI_LIBRARY
)
