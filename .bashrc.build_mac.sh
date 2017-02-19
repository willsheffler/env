function releasebuild {
    pd=`pwd`
    d=$2
    if [ ! $d ]; then d="."; fi
    cd $d
    c1="/usr/bin/c++   -DMAC -DNDEBUG -DCPPDB_EXPORTS -DCPPDB_LIBRARY_PREFIX=\"lib\" -DCPPDB_LIBRARY_SUFFIX=\".dylib\" -DCPPDB_SOVERSION=\"0\" -DCPPDB_DISABLE_THREAD_SAFETY -DCPPDB_DISABLE_SHARED_OBJECT_LOADING -DSQLITE_OMIT_LOAD_EXTENSION -DSQLITE_OMIT_DISABLE_LFS -DSQLITE_THREADSAFE=0 -DCPPDB_WITH_SQLITE3 -pipe -w -ffast-math -funroll-loops -finline-functions -O3 -finline-limit=20000 -s -I$HOME/rosetta_source/cmake/build_release/../.. -I$HOME/rosetta_source/cmake/build_release/../../src -I$HOME/rosetta_source/cmake/build_release/../../external/cxxtest -I$HOME/rosetta_source/cmake/build_release/../../external/boost_1_46_1 -I$HOME/rosetta_source/cmake/build_release/../../external/include -I$HOME/rosetta_source/cmake/build_release/../../src/platform/macos -I$HOME/rosetta_source/cmake/build_release/../../external/dbio -I$HOME/rosetta_source/cmake/build_release/../../external    -pipe -w -o CMakeFiles/$1.dir$HOME/rosetta_source/src/apps/pilot/will/$1.cc.o -c $HOME/rosetta_source/src/apps/pilot/will/$1.cc"
    c2="/usr/bin/c++    -pipe -w -ffast-math -funroll-loops -finline-functions -O3 -finline-limit=20000 -s -Wl,-search_paths_first -Wl,-headerpad_max_install_names   CMakeFiles/$1.dir$HOME/rosetta_source/src/apps/pilot/will/$1.cc.o  -o $1  -L$HOME/rosetta_source/cmake/build_release/../../external/boost_1_46_1 -L$HOME/rosetta_source/cmake/build_release/../../external/lib libdevel.dylib libprotocols.dylib libcore.5.dylib libcore.4.dylib libcore.3.dylib libcore.2.dylib libcore.1.dylib libbasic.dylib libnumeric.dylib libutility.dylib libObjexxFCL.dylib -lz libcppdb.dylib -lz libcppdb.dylib libsqlite3.dylib 
	"
    $c1 && $c2
    if [ $d ]; then cd $pd; fi
}
function releasebuildfast {
    pd=`pwd`
    d=$2
    if [ ! $d ]; then d="."; fi
    cd $d
    c1="/usr/bin/c++   -DMAC -DNDEBUG -DCPPDB_EXPORTS -DCPPDB_LIBRARY_PREFIX=\"lib\" -DCPPDB_LIBRARY_SUFFIX=\".dylib\" -DCPPDB_SOVERSION=\"0\" -DCPPDB_DISABLE_THREAD_SAFETY -DCPPDB_DISABLE_SHARED_OBJECT_LOADING -DSQLITE_OMIT_LOAD_EXTENSION -DSQLITE_OMIT_DISABLE_LFS -DSQLITE_THREADSAFE=0 -DCPPDB_WITH_SQLITE3 -pipe -w -O0 -finline-limit=20000 -s -I$HOME/rosetta_source/cmake/build_release/../.. -I$HOME/rosetta_source/cmake/build_release/../../src -I$HOME/rosetta_source/cmake/build_release/../../external/cxxtest -I$HOME/rosetta_source/cmake/build_release/../../external/boost_1_46_1 -I$HOME/rosetta_source/cmake/build_release/../../external/include -I$HOME/rosetta_source/cmake/build_release/../../src/platform/macos -I$HOME/rosetta_source/cmake/build_release/../../external/dbio -I$HOME/rosetta_source/cmake/build_release/../../external    -pipe -w -o CMakeFiles/$1.dir$HOME/rosetta_source/src/apps/pilot/will/$1.cc.o -c $HOME/rosetta_source/src/apps/pilot/will/$1.cc"
    c2="/usr/bin/c++    -pipe -w -O0 -finline-limit=20000 -s -Wl,-search_paths_first -Wl,-headerpad_max_install_names   CMakeFiles/$1.dir$HOME/rosetta_source/src/apps/pilot/will/$1.cc.o  -o $1  -L$HOME/rosetta_source/cmake/build_release/../../external/boost_1_46_1 -L$HOME/rosetta_source/cmake/build_release/../../external/lib libdevel.dylib libprotocols.dylib libcore.5.dylib libcore.4.dylib libcore.3.dylib libcore.2.dylib libcore.1.dylib libbasic.dylib libnumeric.dylib libutility.dylib libObjexxFCL.dylib -lz libcppdb.dylib -lz libcppdb.dylib libsqlite3.dylib 
	"
    $c1 && $c2
    if [ $d ]; then cd $pd; fi
}

# function releasebuildopt {
#     pd=`pwd`
#     d=$2
#     if [ ! $d ]; then d="."; fi
#     cd $d
#     c1="g++ -o build/src/release/macos/10.7/64/x86/gcc/apps/pilot/will/$1.o -c -isystem external/boost_1_46_1/boost/ -isystem external/boost_1_46_1/boost/ -m64 -march=nocona -mtune=generic -O3 -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s -Wno-unused-variable -DNDEBUG -Isrc -Iexternal/include -Isrc/platform/macos/64/gcc -Isrc/platform/macos/64 -Isrc/platform/macos -Iexternal/boost_1_46_1 -Iexternal/dbio -I/usr/local/include -I/usr/include src/apps/pilot/will/$1.cc"
#     c2="g++ -o build/src/release/macos/10.7/64/x86/gcc/$1.macosgccrelease -m64 -Wl,-stack_size,4000000 build/src/release/macos/10.7/64/x86/gcc/apps/pilot/will/$1.o -Llib -Lexternal/lib -Lbuild/src/release/macos/10.7/64/x86/gcc -Lsrc -Lbuild/external/release/macos/10.7/64/x86/gcc -Lexternal -L/usr/local/lib -L/usr/lib -ldevel -lprotocols -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3"
#     $c1 && $c2
#     if [ $d ]; then cd $pd; fi
# }
# 
# function ompbuild {
#     pd=`pwd`
#     d=$2
#     if [ ! $d ]; then d="."; fi
#     cd $d
#     c1="g++ -o build/src/release/macos/10.7/64/x86/gcc/omp/apps/pilot/will/$1.o -c -isystem external/boost_1_46_1/boost/ -isystem external/boost_1_46_1/boost/ -m64 -march=nocona -mtune=generic -O0 -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s -Wno-unused-variable -fopenmp -DNDEBUG -DUSE_OPENMP -DMULTI_THREADED -Isrc -Iexternal/include -Isrc/platform/macos/64/gcc -Isrc/platform/macos/64 -Isrc/platform/macos -Iexternal/boost_1_46_1 -Iexternal/dbio -I/usr/local/include -I/usr/include src/apps/pilot/will/$1.cc"
#     c2="g++ -o build/src/release/macos/10.7/64/x86/gcc/omp/$1.omp.macosgccrelease -m64 -Wl,-stack_size,4000000 build/src/release/macos/10.7/64/x86/gcc/omp/apps/pilot/will/$1.o -Llib -Lexternal/lib -Lbuild/src/release/macos/10.7/64/x86/gcc/omp -Lsrc -Lbuild/external/release/macos/10.7/64/x86/gcc/omp -Lexternal -L/usr/local/lib -L/usr/lib -ldevel -lprotocols -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3 -lgomp"
#     $c1 && $c2
#     cd $pd
# }
# function ompbuildopt {
#     pd=`pwd`
#     d=$2
#     if [ ! $d ]; then d="."; fi
#     cd $d
#     c1="g++ -o build/src/release/macos/10.7/64/x86/gcc/omp/apps/pilot/will/$1.o -c -isystem external/boost_1_46_1/boost/ -isystem external/boost_1_46_1/boost/ -m64 -march=nocona -mtune=generic -O3 -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s -Wno-unused-variable -fopenmp -DNDEBUG -DUSE_OPENMP -DMULTI_THREADED -Isrc -Iexternal/include -Isrc/platform/macos/64/gcc -Isrc/platform/macos/64 -Isrc/platform/macos -Iexternal/boost_1_46_1 -Iexternal/dbio -I/usr/local/include -I/usr/include src/apps/pilot/will/$1.cc"
#     c2="g++ -o build/src/release/macos/10.7/64/x86/gcc/omp/$1.omp.macosgccrelease -m64 -Wl,-stack_size,4000000 build/src/release/macos/10.7/64/x86/gcc/omp/apps/pilot/will/$1.o -Llib -Lexternal/lib -Lbuild/src/release/macos/10.7/64/x86/gcc/omp -Lsrc -Lbuild/external/release/macos/10.7/64/x86/gcc/omp -Lexternal -L/usr/local/lib -L/usr/lib -ldevel -lprotocols -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3 -lgomp"
#     $c1 && $c2
#     cd $pd
# }
# 
# function debugbuild {
#     pd=`pwd`
#     d=$2
#     if [ ! $d ]; then d="."; fi
#     cd $d
#     c1="g++ -o build/src/debug/macos/10.7/64/x86/gcc/apps/pilot/will/$1.o -c -isystem external/boost_1_46_1/boost/ -isystem external/boost_1_46_1/boost/ -m64 -march=nocona -mtune=generic -O0 -g -ggdb -ffloat-store -Isrc -Iexternal/include -Isrc/platform/macos/64/gcc -Isrc/platform/macos/64 -Isrc/platform/macos -Iexternal/boost_1_46_1 -Iexternal/dbio -I/usr/local/include -I/usr/include src/apps/pilot/will/$1.cc"
#     c2="g++ -o build/src/debug/macos/10.7/64/x86/gcc/$1.macosgccdebug -m64 -Wl,-stack_size,4000000 build/src/debug/macos/10.7/64/x86/gcc/apps/pilot/will/$1.o -Llib -Lexternal/lib -Lbuild/src/debug/macos/10.7/64/x86/gcc -Lsrc -Lbuild/external/debug/macos/10.7/64/x86/gcc -Lexternal -L/usr/local/lib -L/usr/lib -ldevel -lprotocols -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3"
#     $c1 && $c2
#     cd $pd
# }
# 
# function ompbuilddbg {
#     pd=`pwd`
#     d=$2
#     if [ ! $d ]; then d="."; fi
#     cd $d
#     g++ -o build/src/debug/macos/10.7/64/x86/gcc/omp/apps/pilot/will/$1.o -c -isystem external/boost_1_46_1/boost/ -isystem external/boost_1_46_1/boost/ -m64 -march=nocona -mtune=generic -O0 -g -ggdb -ffloat-store -fopenmp -DUSE_OPENMP -DMULTI_THREADED -Isrc -Iexternal/include -Isrc/platform/macos/64/gcc -Isrc/platform/macos/64 -Isrc/platform/macos -Iexternal/boost_1_46_1 -Iexternal/dbio -I/usr/local/include -I/usr/include src/apps/pilot/will/$1.cc
#     g++ -o build/src/debug/macos/10.7/64/x86/gcc/omp/$1.omp.macosgccdebug -m64 -Wl,-stack_size,4000000 build/src/debug/macos/10.7/64/x86/gcc/omp/apps/pilot/will/$1.o -Llib -Lexternal/lib -Lbuild/src/debug/macos/10.7/64/x86/gcc/omp -Lsrc -Lbuild/external/debug/macos/10.7/64/x86/gcc/omp -Lexternal -L/usr/local/lib -L/usr/lib -ldevel -lprotocols -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3 -lgomp
#     cd $pd
# }
# 
