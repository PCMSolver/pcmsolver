#!/usr/bin/env python
import re
import sys

badtests = []
testfail = re.compile(r'^\s*(?P<num>\d+) - (?P<name>\w+(?:-\w+)*) \(Failed\)\s*$')

with open('full_ctest_output.dat', 'r') as outfile:
    ctestout = outfile.readlines()

ctest_exit_status = int(ctestout[0])
if len(ctestout[1:]) == 0:
    sys.stdout.write('\n  <<<  All test cases have passed!  >>>\n\n')
else:
    sys.stdout.write('\n  <<<  Failing outputs follow.  >>>\n\n')

for line in ctestout[1:]:
    linematch = testfail.match(line)
    if linematch:
        bad = linematch.group('name')
        sys.stdout.write('\n\n {} failed. Here is the output:\n'.format(bad))

        badoutfile = 'tests/' + bad + '.log'

        with open(badoutfile, 'r') as ofile:
            sys.stdout.write(ofile.read())

# <<<  return ctest error code  >>>
sys.exit(ctest_exit_status)
