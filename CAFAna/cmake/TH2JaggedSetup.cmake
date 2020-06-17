if(DEFINED USE_TH2JAGGED AND USE_TH2JAGGED)
  include(ExternalProject)

  ExternalProject_Add(TH2Jagged_ext
  PREFIX "${PROJECT_BINARY_DIR}/Ext"
  GIT_REPOSITORY https://github.com/luketpickering/TH2Jagged.git
  GIT_TAG stable
  UPDATE_DISCONNECTED 1
  CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})

  include_directories(${PROJECT_BINARY_DIR}/Ext/src/TH2Jagged_ext/src)

  #Add it after the no-as-needed
  LIST(INSERT ROOT_LIBS 1 TH2Jagged)

  LIST(APPEND EXTRA_LINK_DIRS ${CMAKE_INSTALL_PREFIX}/lib)

  LIST(APPEND EXTRA_CXX_FLAGS -DUSE_TH2JAGGED)

endif()