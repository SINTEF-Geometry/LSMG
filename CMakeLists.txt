PROJECT(lsmg)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)



ADD_SUBDIRECTORY(mba)


# Find modules

FIND_PACKAGE(Boost REQUIRED)


# Include directories

INCLUDE_DIRECTORIES(
  ${lsmg_SOURCE_DIR}/include
  ${mba_SOURCE_DIR}/include
  ${Boost_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  mba
  )


# Make the mba library

FILE(GLOB_RECURSE lsmg_SRCS src/*.C src/*.cpp include/*.h)
ADD_LIBRARY(lsmg ${lsmg_SRCS})


# Set various compiler flags

IF(NOT WIN32)
  # Set warning level to '-Wall' on Linux
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ELSE(NOT WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8
  -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
ENDIF(NOT WIN32)

# Get rid of MESSAGEs in Release mode
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNVERBOSE")


# Apps, examples, tests, ...?

FILE(GLOB lsmg_APPS app/*.C app/*.cpp)
FOREACH(app ${lsmg_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} lsmg ${DEPLIBS})
ENDFOREACH(app)

# Copy data
ADD_CUSTOM_COMMAND(
  TARGET lsmg
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${lsmg_SOURCE_DIR}/Data 
  ${lsmg_BINARY_DIR}/Data
  )


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(lsmg_INSTALL_PREFIX "$ENV{PROGRAMFILES}/LSMG"
    CACHE PATH "Path to install LSMG")
  INSTALL(TARGETS lsmg DESTINATION ${lsmg_INSTALL_PREFIX}/lib)
  INSTALL(DIRECTORY include
    DESTINATION ${lsmg_INSTALL_PREFIX}/include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS lsmg DESTINATION lib/lsmg)
  # include
  INSTALL(DIRECTORY include
    DESTINATION include/lsmg
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)


# CPack stuff
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "LSMG-1.3.0")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
INCLUDE(CPack)
