#!/bin/bash

cmake_build () {
  _BUILD_TOP=$1
  _BUILD_TARGET=$2
  _BUILD_MODE=$3
  
  _CPUS=$(nproc)
  if [ $_CPUS -gt 1 ]
  then
    _MAKE_OPTS="-j $_CPUS"
  fi
  
  export BUILD_DIR=$TOP/out.$_BUILD_TARGET
  mkdir -p $BUILD_DIR
  cd $BUILD_DIR
  case $_BUILD_MODE in
    debug)
      cmake \
        -DCMAKE_TOOLCHAIN_FILE=$_BUILD_TOP/build/$_BUILD_TARGET-toolchain.cmake \
        -DCMAKE_BUILD_TYPE=DEBUG \
        $_BUILD_TOP/build
      ;;
    release)
      cmake \
        -DCMAKE_TOOLCHAIN_FILE=$_BUILD_TOP/build/$_BUILD_TARGET-toolchain.cmake \
        $_BUILD_TOP/build
      ;;
  esac
  return_val=$?
  if [ $return_val = 0 ]; then
    make $_MAKE_OPTS
    return_val=$?
  fi
  cd $_BUILD_TOP
  unset BUILD_DIR
  return $return_val
}

TOP=`pwd`
if [ ! -e $TOP/build ]; then
  echo "run $0 in top directory"
  exit 1
fi

BUILD_TARGET="generic"
BUILD_MODE="debug"
while getopts t:m: OPT
do
  case $OPT in
    t) BUILD_TARGET=$OPTARG
      ;;
    m) BUILD_MODE=$OPTARG
      ;;
  esac
done

case $BUILD_TARGET in
  generic)
    cmake_build $TOP $BUILD_TARGET $BUILD_MODE
    ;;
  stm32*)
    ./build.$BUILD_TARGET.sh
    ;;
esac
