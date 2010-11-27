#!/usr/bin/python
import sys,os,random

list = open(sys.argv[1]).readlines()
random.shuffle(list)
print "".join(list)
