
'''
Script checks whether any file under library_path is newer than timestamp_file.
If this is the case, timestamp_path is removed.
'''

import os
import sys
import fnmatch
import shutil

timestamp_file = sys.argv[1]
timestamp_path = sys.argv[2]
library_path   = sys.argv[3]

# silently exit if time stamp file does not exist
if not os.path.isfile(timestamp_file):
    sys.exit(0)

# get time of configure stamp
stamp_time = os.path.getmtime(timestamp_file)

# find out whether any file within library_path is newer than timestamp_file
reset_stamp = False
for root, dirnames, filenames in os.walk(library_path):
    if reset_stamp:
        break
    for filename in fnmatch.filter(filenames, '*'):
        f = os.path.join(root, filename)
        if os.path.getmtime(f) > stamp_time:
            reset_stamp = True
            break

if reset_stamp:
    # sanity check
    if os.path.dirname(timestamp_file) == timestamp_path:
        shutil.rmtree(timestamp_path)
