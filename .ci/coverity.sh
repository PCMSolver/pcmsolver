. "$CI_SCRIPTS/common.sh"

# temporarily disable error checking, the coverity script exits with
# status code 1 whenever it (1) fails OR (2) is not on the correct
# branch.
set +e
curl https://scan.coverity.com/scripts/travisci_build_coverity_scan.sh 
sed -i -e 's@https:\/\/scan.coverity.com\/download\/@https:\/\/scan.coverity.com\/download\/cxx\/@g' travisci_build_coverity_scan.sh  
COVERITY_SCAN_PROJECT_NAME="PCMSolver/pcmsolver" \
COVERITY_SCAN_NOTIFICATION_EMAIL="roberto.d.remigio@uit.no" \
COVERITY_SCAN_BRANCH_PATTERN="release" \
COVERITY_SCAN_BUILD_COMMAND_PREPEND="python setup --cxx=g++ --cc=gcc --fc=gfortran --type=debug --tests --coverage" \
COVERITY_SCAN_BUILD_COMMAND="make" \
bash travisci_build_coverity_scan.sh
set -e

exit 0
