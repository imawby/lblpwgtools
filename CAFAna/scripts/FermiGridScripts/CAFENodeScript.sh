#!/bin/bash

LOG_TO_IFDH=0

LOGYLOG () {
  echo ${1}
  if [ ${LOG_TO_IFDH} == "1" ]; then
    ifdh log ${1}
  fi
}

if [ ! -z ${2} ]; then
  LOG_TO_IFDH=${2}
  if [ ${LOG_TO_IFDH} == "1" ]; then
    LOGYLOG "[INFO]: Logging to ifdh log."
  fi
fi

#if [ -z ${INPUT_TAR_FILE} ]; then
#  LOGYLOG "[ERROR]: Expected to recieve an input file."
#  exit 1
#fi

if [ ! -e ${INPUT_TAR_DIR_LOCAL}/CAFAna/CAFECommands.cmd ]; then
  LOGYLOG "[ERROR]: Expected to recieve a command file @ CAFAna/CAFECommands.cmd but didn't."
  ls ${INPUT_TAR_DIR_LOCAL}/CAFAna
  exit 2
fi

PNFS_PATH_APPEND=${1}
if [ -z ${1} ]; then
  LOGYLOG "[ERROR]: Failed to find PNFS_PATH_APPEND passed on command line."
  exit 3
fi

printenv

set -x #start bash debugging at this point
LOGYLOG "Start $(date)"
LOGYLOG "Site:${GLIDEIN_ResourceName}"
LOGYLOG "The worker node is " `hostname` "OS: " `uname -a`
LOGYLOG "The user id is $(whoami)"
LOGYLOG "The output of id is: $(id)"
set +x #stop bash debugging at this point

if [ -z ${GRID_USER} ]; then
  GRID_USER=$(basename $X509_USER_PROXY | cut -d "_" -f 2)
fi

if [ -z ${GRID_USER} ]; then
  LOGYLOG "Failed to get GRID_USER."
  exit 4
fi

LOGYLOG "GRID_USER is ${GRID_USER}"

mv ${INPUT_TAR_DIR_LOCAL}/CAFAna $_CONDOR_SCRATCH_DIR/

cd $_CONDOR_SCRATCH_DIR

export CAFANA=${INPUT_TAR_DIR_LOCAL}/CAFAna
#source ${CAFANA}/CAFAnaEnv.sh
source ${CAFANA}/CAFAnaEnv.sh
LOGYLOG "BOOST_DIR IS ${BOOST_DIR}"
LOGYLOG "BOOST_INC IS ${BOOST_INC}"
LOGYLOG "BOOST_LIB IS ${BOOST_LIB}"
LOGYLOG "ROOT_INCLUDE_PATH IS ${ROOT_INCLUDE_PATH}"
LOGYLOG "LD_LIBRARY_PATH IS ${LD_LIBRARY_PATH}"
LOGYLOG "SETUP_BOOST IS ${SETUP_BOOST}"
LOGYLOG "BOOST_LIB IS {$BOOST_LIB}"
LOGYLOG "BOOST_FQ_DIR IS ${BOOST_FQ_DIR}"
LOGYLOG "PKG_CONFIG_PATH IS ${PKG_CONFIG_PATH}"
LOGYLOG "CMAKE_PREFIX_PATH IS ${CMAKE_PREFIX_PATH}"

voms-proxy-info --all
setup ifdhc
ups active
kx509

echo "X509_USER_PROXY is"
echo $X509_USER_PROXY

#setup_fnal_security

export IFDH_CP_UNLINK_ON_ERROR=1;
export IFDH_CP_MAXRETRIES=2;

PNFS_OUTDIR_STUB=/pnfs/dune/persistent/users/${GRID_USER}/${PNFS_PATH_APPEND}
LOGYLOG "Output stub dir is ${PNFS_OUTDIR_STUB}"

ifdh ls ${PNFS_OUTDIR_STUB}

if [ $? -ne 0 ]; then
  LOGYLOG "Unable to read ${PNFS_OUTDIR_STUB}. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
  exit 5
fi

PNFS_OUTDIR=${PNFS_OUTDIR_STUB}/${CLUSTER}.${PROCESS}/

ifdh mkdir ${PNFS_OUTDIR}
ifdh ls ${PNFS_OUTDIR}
if [ $? -ne 0 ]; then
  LOGYLOG "Unable to make ${PNFS_OUTDIR}."
  exit 2
fi

LOGYLOG "Output dir is ${PNFS_OUTDIR}"

(( LINE_N = ${PROCESS} + 1 ))

LINE=$(cat ${CAFANA}/CAFECommands.cmd | head -${LINE_N} | tail -1)

LOGYLOG "Running command: ${LINE}"

SCRIPT_NAME=$(echo ${LINE} | cut -f 1 -d " ")
OUTPUT_NAME=$(echo ${LINE} | cut -f 3 -d " ")
N_THROWS=$(echo ${LINE} | cut -f 4 -d " ")
N_TEST_POINTS=$(echo ${LINE} | cut -f 5 -d " ")
TEST_VALUE_INDEX=$(echo ${LINE} | cut -f 6 -d " ")

LOGYLOG "Trying to open tarball"

ls ${CAFANA}

LOGYLOG "Maybe opened tarball? who knows?"

LOGYLOG "Running script ${SCRIPT_NAME} and expecting output ${OUTPUT_NAME}"

if [ ! -e  ${CAFANA}/scripts/${SCRIPT_NAME} ]; then
  LOGYLOG "[ERROR]: Failed to find expected script: ${CAFANA}/scripts/${SCRIPT_NAME}"
  exit 6
fi

#cp ${CAFANA}/scripts/common_fit_definitions.C .
#cp ${CAFANA}/scripts/${SCRIPT_NAME} .

LOGYLOG "Running script @ $(date)"

#LOGYLOG "cafe -q -b ${CAFANA}/scripts/${SCRIPT_NAME} $(echo ${LINE} | cut -f 2- -d " ") &> ${SCRIPT_NAME}.log"
#cafe -q -b ${CAFANA}/scripts/${SCRIPT_NAME} $(echo ${LINE} | cut -f 2- -d " ") &> ${SCRIPT_NAME}.log

ups active

unsetup python
setup python v3_7_2
setup gcc v9_3_0
setup boost v1_75_0  
setup dune_oslibs
unsetup ifdh
setup ifdhc v2_5_16 -q e20:p392:prof

LOGYLOG "IFDHC_LIB IS ${IFDHC_LIB}"
LOGYLOG "BOOST_DIR IS ${BOOST_DIR}"
LOGYLOG "BOOST_INC IS ${BOOST_INC}"
LOGYLOG "BOOST_LIB IS ${BOOST_LIB}"

FILE_NAME=${CAFANA}/StateFilesAllSystematicsSplitBySign.root

chmod ugo+rwx ${FILE_NAME}

ls -lth ${CAFANA}
ls -lth ${CAFANA}/scripts

pwd
cafe -q -b ${CAFANA}/scripts/${SCRIPT_NAME} ${N_THROWS} ${N_TEST_POINTS} ${TEST_VALUE_INDEX} &> ${SCRIPT_NAME}.log  

setup ifdhc

LOGYLOG "cafe ${CAFANA}/scripts/${SCRIPT_NAME} &> ${SCRIPT_NAME}.log"

LOGYLOG "Copying output @ $(date)"

LOGYLOG "ifdh cp -D $IFDH_OPTION ${SCRIPT_NAME}.log ${PNFS_OUTDIR}/"

ifdh cp -D $IFDH_OPTION ${SCRIPT_NAME}.log ${PNFS_OUTDIR}/

if [ ! -e ${OUTPUT_NAME} ]; then
  LOGYLOG "[WARN]: Failed to produce expected output file (${OUTPUT_NAME})."
fi

for f in *.root; do
  LOGYLOG "ifdh cp -D $IFDH_OPTION ${f} ${PNFS_OUTDIR}/"
  ifdh cp -D $IFDH_OPTION ${f} ${PNFS_OUTDIR}/
done

echo "All stop @ $(date)"
ifdh log "All stop @ $(date)"
