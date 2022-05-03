#!/bin/bash

#CMAKE_BUILD_TYPE="Release"
#CMAKE_BUILD_TYPE="Debug"

DIR="./build/"
if [ -d "$DIR" ]; then
  # Take action if $DIR exists. #
  echo "building project in ${DIR}..."
else
  mkdir build
  echo "mkdir ${DIR} ... "
fi

echo "build dir: $DIR"
g++ -std=c++17 -O3 -finline-functions Isotropic_emission.cpp -o build/isotropic.exe;
#cmake . -B${DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}; cd ${DIR}; make VERBOSE=1


