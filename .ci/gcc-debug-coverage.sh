. "$CI_SCRIPTS/common.sh"

pip install cpp-coveralls --user `whoami`

python setup.py --cxx=g++ --cc=gcc --fc=gfortran --type=debug --tests --coverage
cd build
make
ctest

coveralls --root .. -E ".*external.*" -E ".*CMakeFiles.*" -E ".*generated.*" -E ".*tests*" -E ".*tools.*" -E ".*cmake*" -E ".*doc*" -E ".*examples*" || echo 'coveralls upload failed.'
