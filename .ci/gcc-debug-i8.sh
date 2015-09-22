. "$CI_SCRIPTS/common.sh"

python setup.py --cxx=g++ --cc=gcc --fc=gfortran --type=debug --tests --int64
cd build
make
ctest
