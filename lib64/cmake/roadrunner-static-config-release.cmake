#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "roadrunner-static::roadrunner-static" for configuration "Release"
set_property(TARGET roadrunner-static::roadrunner-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(roadrunner-static::roadrunner-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libroadrunner-static.a"
  )

list(APPEND _cmake_import_check_targets roadrunner-static::roadrunner-static )
list(APPEND _cmake_import_check_files_for_roadrunner-static::roadrunner-static "${_IMPORT_PREFIX}/lib64/libroadrunner-static.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
