SET(TARGET test15_cad_bathymetry)
SET(_main test15_cad_bathymetry.cc)


#set(OPENCASCADE_DIR /home/ole/dev/tools/oce-0.18.3)
set(OPENCASCADE_DIR /usr/local/oce-upstream-V7_5_0)
message("OpenCASCADE directory:")
message(${OPENCASCADE_DIR})

# Set the _main variable to empty if you don't want an executable
# but only a library
#SET(_main "")

############################################################
# Normally you shouldn't need to change anything below.
############################################################
# Declare all source files the target consists of:
file(GLOB _files source/*cc)

# Don't compile the main file into the library
IF(NOT "${_main}" STREQUAL "")
    LIST(REMOVE_ITEM _files
      ${CMAKE_CURRENT_SOURCE_DIR}/${_main}
    )
ENDIF()


CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

FIND_PACKAGE(deal.II 8.3 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)

IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING 
       "Choose the type of build, options are: Debug, Release or DebugRelease"
       )
ELSE()
   SET(CMAKE_BUILD_TYPE ${DEAL_II_BUILD_TYPE} CACHE STRING
       "Choose the type of build, options are: Debug, Release or DebugRelease"
       )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()


PROJECT(${TARGET})

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
set(CMAKE_MACOSX_RPATH 1)
endif()


FIND_PACKAGE(deal2lkit 1.0 REQUIRED
  HINTS ${D2K_DIR} $ENV{D2K_DIR} $ENV{DEAL2LKIT_DIR}
  )
D2K_INITIALIZE_CACHED_VARIABLES()

FIND_PACKAGE(Boost 1.56 REQUIRED COMPONENTS filesystem program_options)

# We one library and one target for each type of deal.II library
# we found. If you compiled deal.II with both Release and Debug
# mode, this will generate both Release and Debug programs for you
# The debug library and program are postfixed with ".g"
SET(_d2_build_types "Release;Debug")
SET(Release_postfix "")
SET(Debug_postfix ".g")

FOREACH(_build_type ${_d2_build_types})
    # Postfix to use everywhere
    SET(_p "${${_build_type}_postfix}")
    # Only build this type, if deal.II was compiled with it.
    IF(CMAKE_BUILD_TYPE MATCHES "${_build_type}" AND
       DEAL_II_BUILD_TYPE MATCHES "${_build_type}"  AND
       D2K_BUILD_TYPE MATCHES "${_build_type}")
	MESSAGE("-- Found ${_build_type} version of deal.II.")
	MESSAGE("-- Found ${_build_type} version of deal2lkit.")


	STRING(TOUPPER "${_build_type}" _BUILD_TYPE)

     FOREACH(_dim 3)
		SET(_exe "${TARGET}_${_dim}d${${_build_type}_postfix}")
		MESSAGE("-- Configuring executable ${_exe}")
		ADD_EXECUTABLE(${_exe} ${_main})
		target_compile_definitions(${_exe} PUBLIC DEAL_II_DIMENSION=${_dim})
		TARGET_LINK_LIBRARIES(${_exe} 
              ${_lib} 
              Boost::program_options Boost::filesystem
${OPENCASCADE_DIR}/lib/libTKernel.so
${OPENCASCADE_DIR}/lib/libTKMath.so
${OPENCASCADE_DIR}/lib/libTKG2d.so
${OPENCASCADE_DIR}/lib/libTKG3d.so
${OPENCASCADE_DIR}/lib/libTKGeomBase.so
${OPENCASCADE_DIR}/lib/libTKBRep.so
${OPENCASCADE_DIR}/lib/libTKGeomAlgo.so
${OPENCASCADE_DIR}/lib/libTKTopAlgo.so
${OPENCASCADE_DIR}/lib/libTKPrim.so
${OPENCASCADE_DIR}/lib/libTKBO.so
${OPENCASCADE_DIR}/lib/libTKShHealing.so
${OPENCASCADE_DIR}/lib/libTKBool.so
${OPENCASCADE_DIR}/lib/libTKHLR.so
${OPENCASCADE_DIR}/lib/libTKFillet.so
${OPENCASCADE_DIR}/lib/libTKOffset.so
${OPENCASCADE_DIR}/lib/libTKFeat.so
${OPENCASCADE_DIR}/lib/libTKMesh.so
${OPENCASCADE_DIR}/lib/libTKService.so
${OPENCASCADE_DIR}/lib/libTKV3d.so
${OPENCASCADE_DIR}/lib/libTKOpenGl.so
${OPENCASCADE_DIR}/lib/libTKMeshVS.so
${OPENCASCADE_DIR}/lib/libTKCDF.so
${OPENCASCADE_DIR}/lib/libTKLCAF.so
${OPENCASCADE_DIR}/lib/libTKCAF.so
${OPENCASCADE_DIR}/lib/libTKBinL.so
${OPENCASCADE_DIR}/lib/libTKXmlL.so
${OPENCASCADE_DIR}/lib/libTKBin.so
${OPENCASCADE_DIR}/lib/libTKXml.so
${OPENCASCADE_DIR}/lib/libTKTObj.so
${OPENCASCADE_DIR}/lib/libTKBinTObj.so
${OPENCASCADE_DIR}/lib/libTKXmlTObj.so
${OPENCASCADE_DIR}/lib/libTKXSBase.so
${OPENCASCADE_DIR}/lib/libTKSTEPBase.so
${OPENCASCADE_DIR}/lib/libTKSTEPAttr.so
${OPENCASCADE_DIR}/lib/libTKSTEP209.so
${OPENCASCADE_DIR}/lib/libTKSTEP.so
${OPENCASCADE_DIR}/lib/libTKIGES.so
${OPENCASCADE_DIR}/lib/libTKXDEIGES.so
${OPENCASCADE_DIR}/lib/libTKXCAF.so
${OPENCASCADE_DIR}/lib/libTKXDESTEP.so
${OPENCASCADE_DIR}/lib/libTKSTL.so
${OPENCASCADE_DIR}/lib/libTKXmlXCAF.so
${OPENCASCADE_DIR}/lib/libTKVRML.so
${OPENCASCADE_DIR}/lib/libTKBinXCAF.so
              )
		D2K_SETUP_TARGET(${_exe} ${_BUILD_TYPE})
 	ENDFOREACH()

  ENDIF()
ENDFOREACH()

