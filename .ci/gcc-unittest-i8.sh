. "$CI_SCRIPTS/common.sh"

sudo pip install cpp-coveralls

python setup --cxx=g++ --cc=gcc --fc=gfortran --type=debug --tests --int64
cd build
make
make test
