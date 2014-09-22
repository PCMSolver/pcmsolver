. "$CI_SCRIPTS/common.sh"

python setup --cxx=clang++ --cc=clang --fc=gfortran --type=debug --tests
cd build
$MAKE_CMD
$MAKE_CMD test
