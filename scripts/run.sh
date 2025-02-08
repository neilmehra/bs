#/bin/sh
if [ ! -d ./build ]; then
  echo "Build dir doesn't exist."
  scripts/build.sh
else
  cd build
  cmake --build .
  cd ..
  scripts/test.sh
fi

