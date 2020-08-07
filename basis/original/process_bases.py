#!/usr/bin/python

import glob
from subprocess import call

for name in glob.glob('*.dat'):
	call(["./process_file",name])
