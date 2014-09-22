. "$CI_SCRIPTS/common.sh"

sudo pip install cpp-coveralls

python setup --cxx=g++ --cc=gcc --fc=gfortran --type=debug --tests --coverage
cd build
make
make test

coveralls --root .. -E ".*external.*" -E ".*CMakeFiles.*" -E ".*generated.*" -E ".*tests*" -E ".*tools.*" || echo 'coveralls upload failed.'
