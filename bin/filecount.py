#!/usr/bin/env python
import sys,os
for dir_in in sys.argv[1:]:
	fc = {}
	for l in os.popen("find %s"%dir_in).readlines():
		l = l.replace("//","/")
		paths = l.split("/")
		path = paths[0]
		if not path in fc: fc[path] = 0
		fc[path] += 1
		for n in paths[1:-1]:
			if not n: continue
			path += "/"+n
			if not path in fc: fc[path] = 0
			fc[path] += 1

dirs = fc.keys()
dirs.sort()
for d in dirs:
	print "%7i %s"%(fc[d],d)

		