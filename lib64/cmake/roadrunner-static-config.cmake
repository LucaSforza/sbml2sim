# Generated by CMake

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.8)
   message(FATAL_ERROR "CMake >= 2.8.3 required")
endif()
if(CMAKE_VERSION VERSION_LESS "2.8.3")
   message(FATAL_ERROR "CMake >= 2.8.3 required")
endif()
cmake_policy(PUSH)
cmake_policy(VERSION 2.8.3...3.29)
#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Protect against multiple inclusion, which would fail when already imported targets are added once more.
set(_cmake_targets_defined "")
set(_cmake_targets_not_defined "")
set(_cmake_expected_targets "")
foreach(_cmake_expected_target IN ITEMS roadrunner-static::roadrunner-static)
  list(APPEND _cmake_expected_targets "${_cmake_expected_target}")
  if(TARGET "${_cmake_expected_target}")
    list(APPEND _cmake_targets_defined "${_cmake_expected_target}")
  else()
    list(APPEND _cmake_targets_not_defined "${_cmake_expected_target}")
  endif()
endforeach()
unset(_cmake_expected_target)
if(_cmake_targets_defined STREQUAL _cmake_expected_targets)
  unset(_cmake_targets_defined)
  unset(_cmake_targets_not_defined)
  unset(_cmake_expected_targets)
  unset(CMAKE_IMPORT_FILE_VERSION)
  cmake_policy(POP)
  return()
endif()
if(NOT _cmake_targets_defined STREQUAL "")
  string(REPLACE ";" ", " _cmake_targets_defined_text "${_cmake_targets_defined}")
  string(REPLACE ";" ", " _cmake_targets_not_defined_text "${_cmake_targets_not_defined}")
  message(FATAL_ERROR "Some (but not all) targets in this export set were already defined.\nTargets Defined: ${_cmake_targets_defined_text}\nTargets not yet defined: ${_cmake_targets_not_defined_text}\n")
endif()
unset(_cmake_targets_defined)
unset(_cmake_targets_not_defined)
unset(_cmake_expected_targets)


# Compute the installation prefix relative to this file.
get_filename_component(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
if(_IMPORT_PREFIX STREQUAL "/")
  set(_IMPORT_PREFIX "")
endif()

# Create imported target roadrunner-static::roadrunner-static
add_library(roadrunner-static::roadrunner-static STATIC IMPORTED)

set_target_properties(roadrunner-static::roadrunner-static PROPERTIES
  INTERFACE_COMPILE_DEFINITIONS "STATIC_RR"
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include"
  INTERFACE_LINK_LIBRARIES "pthread;dl;m;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMOrcJIT.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMPasses.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMObjCARCOpts.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMCoroutines.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMipo.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMInstrumentation.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMVectorize.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMLinker.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMIRReader.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMAsmParser.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMFrontendOpenMP.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMJITLink.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMX86Disassembler.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMX86AsmParser.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMX86CodeGen.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMCFGuard.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMGlobalISel.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMX86Desc.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMX86Info.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMMCDisassembler.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMSelectionDAG.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMAsmPrinter.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMDebugInfoMSF.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMDebugInfoDWARF.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMCodeGen.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMScalarOpts.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMInstCombine.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMAggressiveInstCombine.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMTransformUtils.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMBitWriter.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMMCJIT.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMExecutionEngine.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMTarget.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMAnalysis.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMProfileData.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMRuntimeDyld.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMOrcTargetProcess.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMOrcShared.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMObject.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMTextAPI.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMMCParser.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMBitReader.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMMC.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMDebugInfoCodeView.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMCore.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMRemarks.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMBitstreamReader.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMBinaryFormat.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMSupport.a;/__w/roadrunner/llvm-13.x-binaries/lib/libLLVMDemangle.a;rr-libstruct::rr-libstruct-static;sbml-static;zlib::zlibstatic;nleq2::nleq2-static;nleq1::nleq1-static;clapack::lapack;clapack::blas;clapack::f2c;SUNDIALS::generic_static;SUNDIALS::nvecserial_static;SUNDIALS::nvecmanyvector_static;SUNDIALS::sunmatrixband_static;SUNDIALS::sunmatrixdense_static;SUNDIALS::sunmatrixsparse_static;SUNDIALS::sunlinsolband_static;SUNDIALS::sunlinsoldense_static;SUNDIALS::sunlinsolpcg_static;SUNDIALS::sunlinsolspbcgs_static;SUNDIALS::sunlinsolspfgmr_static;SUNDIALS::sunlinsolspgmr_static;SUNDIALS::sunlinsolsptfqmr_static;SUNDIALS::sunnonlinsolnewton_static;SUNDIALS::sunnonlinsolfixedpoint_static;SUNDIALS::cvodes_static;SUNDIALS::kinsol_static;Poco::Foundation;Poco::Net;Poco::XML;expat::expat"
)

# Load information for each installed configuration.
file(GLOB _cmake_config_files "${CMAKE_CURRENT_LIST_DIR}/roadrunner-static-config-*.cmake")
foreach(_cmake_config_file IN LISTS _cmake_config_files)
  include("${_cmake_config_file}")
endforeach()
unset(_cmake_config_file)
unset(_cmake_config_files)

# Cleanup temporary variables.
set(_IMPORT_PREFIX)

# Loop over all imported files and verify that they actually exist
foreach(_cmake_target IN LISTS _cmake_import_check_targets)
  if(CMAKE_VERSION VERSION_LESS "3.28"
      OR NOT DEFINED _cmake_import_check_xcframework_for_${_cmake_target}
      OR NOT IS_DIRECTORY "${_cmake_import_check_xcframework_for_${_cmake_target}}")
    foreach(_cmake_file IN LISTS "_cmake_import_check_files_for_${_cmake_target}")
      if(NOT EXISTS "${_cmake_file}")
        message(FATAL_ERROR "The imported target \"${_cmake_target}\" references the file
   \"${_cmake_file}\"
but this file does not exist.  Possible reasons include:
* The file was deleted, renamed, or moved to another location.
* An install or uninstall procedure did not complete successfully.
* The installation package was faulty and contained
   \"${CMAKE_CURRENT_LIST_FILE}\"
but not all the files it references.
")
      endif()
    endforeach()
  endif()
  unset(_cmake_file)
  unset("_cmake_import_check_files_for_${_cmake_target}")
endforeach()
unset(_cmake_target)
unset(_cmake_import_check_targets)

# This file does not depend on other imported targets which have
# been exported from the same project but in a separate export set.

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
cmake_policy(POP)
