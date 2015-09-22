. "$CI_SCRIPTS/common.sh"

python setup.py --cxx=clang++ --cc=clang --fc=gfortran --type=debug --tests
cd build
make
ctest
