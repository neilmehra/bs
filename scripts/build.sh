#/bin/sh
export CC=/usr/bin/clang
export CXX=/usr/bin/clang++
rm -rf build
mkdir build
cd build
cmake -DTBB_TEST=OFF ..
cmake --build .
