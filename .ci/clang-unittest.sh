. "$CI_SCRIPTS/common.sh"

python setup --cxx=clang++ --cc=clang --fc=gfortran --type=debug --tests
cd build
make
make test
