#!/usr/bin/env python
import sys,os,random,glob

nfiles = int(sys.argv[1])

files = sys.argv[2:]
if len(files) is 1:
	print "globbing",files[0]
	files = glob.glob(files[0])

random.shuffle(files)
if len(files) <= nfiles:
	print "doing nothing, too few files to prune to requested number"

if not os.path.exists("__REMOVED_FILES__"): os.mkdir("__REMOVED_FILES__")
for f in files[nfiles:]:
	newpath = "__REMOVED_FILES__/"+os.path.dirname(f)
	newname = "__REMOVED_FILES__/"+f
	if not os.path.exists(newpath):
		os.mkdir(newpath)
	print "moving file",f,"to",newname
	os.rename(f,newname)
