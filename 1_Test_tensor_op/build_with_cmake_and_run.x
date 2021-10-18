#!/bin/bash -f

if [[ -d build ]]
then
    rm -rf build
fi

mkdir -p build
cd build
cmake .. -DCMAKE_PREFIX_PATH=$PFUNIT_DIR 
make
./my_test_tensor_op1
#./my_test_tensor_llpack

# or
# ctest --verbose

