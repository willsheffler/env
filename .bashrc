#!/bin/bash

source ~/.alias.bash # figure out how to do this right
source ~/.alias      # figure out how to do this right
source ~/.minirosetta

export EDITOR="mate -w"
export SVN="https://svn.rosettacommons.org/source/"
export PATH=/opt/local/bin:/Applications:~/bin:/usr/local/bin:/usr/local/sbin:/usr/local/mysql/bin:/opt/local/sbin:$PATH
export PS1="\h:\W/> " 

function buildappdbg {
	cmd="g++ -o build/src/debug/macos/10.6/64/x86/gcc/apps/pilot/will/$1.o -c -pipe -ffor-scope -W -Wall -pedantic -Wno-long-long -m64 -march=nocona -mtune=generic -O0 -g -ggdb -ffloat-store -Isrc -Iexternal/include -Isrc/platform/macos/64/gcc -Isrc/platform/macos/64 -Isrc/platform/macos -Iexternal/boost_1_38_0 -I/usr/local/include -I/usr/include src/apps/pilot/will/$1.cc"
	echo $cmd
	$cmd
	cmd="g++ -o build/src/debug/macos/10.6/64/x86/gcc/$1.macosgccdebug -m64 -Wl,-stack_size,4000000 build/src/debug/macos/10.6/64/x86/gcc/apps/pilot/will/$1.o -Llib -Lexternal/lib -Lbuild/src/debug/macos/10.6/64/x86/gcc -Lsrc -L/usr/lib -ldevel -lprotocols -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz"
	echo $cmd
	$cmd
}

function buildapp {
	cmd="g++ -o build/src/release/macos/10.6/64/x86/gcc/apps/pilot/will/$1.o -c -pipe -ffor-scope -W -Wall -pedantic -Wno-long-long -m64 -march=nocona -mtune=generic -O3 -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s -Wno-unused-variable -DNDEBUG -Isrc -Iexternal/include -Isrc/platform/macos/64/gcc -Isrc/platform/macos/64 -Isrc/platform/macos -Iexternal/boost_1_38_0 -I/usr/local/include -I/usr/include src/apps/pilot/will/$1.cc"
	echo $cmd
	$cmd
	cmd="g++ -o build/src/release/macos/10.6/64/x86/gcc/$1.macosgccrelease -m64 -Wl,-stack_size,4000000 build/src/release/macos/10.6/64/x86/gcc/apps/pilot/will/$1.o -Llib -Lexternal/lib -Lbuild/src/release/macos/10.6/64/x86/gcc -Lsrc -L/usr/lib -ldevel -lprotocols -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz"
	echo $cmd
	$cmd
}

