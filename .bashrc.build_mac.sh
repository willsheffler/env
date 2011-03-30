function buildapp {
	pd=`pwd`
	d=$2
	if [ ! $d ]; then d="."; fi
	cd $d
	g++ -o build/src/release/macos/10.6/64/x86/gcc/apps/pilot/will/$1.o -c -pipe -ffor-scope -W -Wall -pedantic -Wno-long-long -m64 -march=nocona -mtune=generic -O3 -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s -Wno-unused-variable -DNDEBUG -Isrc -Iexternal/include -Isrc/platform/macos/64/gcc -Isrc/platform/macos/64 -Isrc/platform/macos -Iexternal/boost_1_38_0 -I/usr/local/include -I/usr/include src/apps/pilot/will/$1.cc
	g++ -o build/src/release/macos/10.6/64/x86/gcc/$1.macosgccrelease -m64 -Wl,-stack_size,4000000 build/src/release/macos/10.6/64/x86/gcc/apps/pilot/will/$1.o -Lexternal/lib -Lbuild/src/release/macos/10.6/64/x86/gcc -Lsrc -L/usr/lib -ldevel -lprotocols -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz
	cd $pd
}