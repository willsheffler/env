import string
from glob import glob
from os import popen,system,chdir,remove,getcwd
from os.path import exists
from popen2 import popen2
from math import floor,sqrt
from operator import add
from sys import stderr,argv,exit
from whrandom import random

def run(command,safe=0):
    print command
    if not safe:
        system(command)

def log(s):
    stderr.write(s)
    if s[-1] != '\n':
        stderr.write('\n')

def mkdir(dir):
    if not exists(dir):
        system('mkdir '+dir)

def Increment( D, key, count=1):
    if not D.has_key(key):
        D[key] = 0
    D[key] = D[key] + count
