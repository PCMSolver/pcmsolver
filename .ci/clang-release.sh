. "$CI_SCRIPTS/common.sh"

python setup --cxx=clang++ --cc=clang --fc=gfortran --type=release --tests
cd build
make
ctest
