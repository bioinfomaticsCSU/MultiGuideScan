This directory contains example scripts to test run the GuideScan software.

/batch and /chromFa are data directories that the two bash scripts in this directory access in order to run the GuideScan software.

guidescan-example-yeast-practical.sh is a bash script that will execute several use cases of the GuideScan software. There are two commands that are presently commented out in this script. These commands fall under the titles:

"Insert Rule Set 2 cutting efficiency scores"
"Transfer Rule Set 2 cutting efficiency scores from original database to more detailed database"

These commands can be uncommented if the user's system has the python module sklearn version 0.16.1 which can be found here: https://pypi.python.org/pypi/scikit-learn/0.16.1

guidescan-example-ebola-fast.sh is the exact same as uniguide-example-yeast-practical.sh except it processes the ebola genome. This script is much faster and its purpose is to serve as an integrative test for the GuideScan software
