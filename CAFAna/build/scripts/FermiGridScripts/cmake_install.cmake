# Install script for directory: /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts/FermiGridScripts

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/Linux")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "DEBUG")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts/FermiGridScripts" TYPE PROGRAM FILES
    "/dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts/FermiGridScripts/tarball.sh"
    "/dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts/FermiGridScripts/FarmBuildInterps.sh"
    "/dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts/FermiGridScripts/BuildInterps.sh"
    "/dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts/FermiGridScripts/FarmCAFENodeScript.sh"
    "/dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/scripts/FermiGridScripts/CAFENodeScript.sh"
    )
endif()

