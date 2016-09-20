import os
import sys
from utils import *

newlist = map_identifiers(sys.argv[1], sys.argv[2], '\t', '\t', 0, 1, True)
writeLIST_to_file(newlist, sys.argv[3])

