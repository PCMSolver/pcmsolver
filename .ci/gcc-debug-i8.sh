. "$CI_SCRIPTS/common.sh"

python setup.py --cxx=g++ --cc=gcc --fc=gfortran --type=debug --int64
cd build
make
ctest
