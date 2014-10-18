. "$CI_SCRIPTS/common.sh"

python setup --cxx=clang++ --cc=clang --fc=gfortran --type=release --tests --disable_python_embedding
cd build
make
make test
