function setup_relax_native {
    source=$1
    dest=$2
    i=$3
    j=$4
    j=${j:0:4}
    echo "making relax native job setup for:" $i $j
    #if [ -e $source/$i/$j/rosetta_data/*0001.pdb ]; then
    #    echo "    `/bin/ls $source/$i/$j/rosetta_data/*0001.pdb`"
    #else 
    echo "   finding native"
    native="$source/$i/$j/rosetta_data/$j.pdb"
    if [ ! -e $native ]; then
	#native="`/bin/ls $source/$i/$j/rosetta_data/$j..pdb`"
	#if [ ! -e $native ]; then
	    native=""
	    #fi
    fi
    if [ -z $native ]; then
	echo "WARNING: no native found for:" $i $j
	continue
    fi
    echo "    generating idealized native"
    tmpdir=`pwd`
    cd $source/$i/$j/rosetta_data
    nice ~/rosetta.gcc -idealize -s $native -paths ~/paths/paths.absolute.txt > idealize.log 2> idealize.err
    if [ ! -e ${j}_0001.pdb ]; then
	echo "WARNING: no idealized native for:" $i $j
	continue
    fi
    echo "    copying rosetta data"
    mkdir -p $dest/$i/$j/
    mkdir -p $dest/$i/$j/relax_native
    cp ~/paths/paths.relaxnative.txt $dest/$i/$j/paths.txt
	nice cp -r $source/$i/$j/rosetta_data $dest/$i/$j/rosetta_data
    cd $dest/$i/$j
    chain=`/bin/ls rosetta_data | egrep "..$j.03......._v..3"`
    chain=${chain:6:1}
    #echo  $j $chain
    echo "    making condor.input:" rn $j $chain "-s" ${j}_0001.pdb
    ~/scripts/condor/make_relax_condor_input.py rn $j $chain rosetta_data/${j}_0001.pdb > condor.input 2> make_relax_condor_input.err
    cd $tmpdir
}

listfile="$1"
listerr="$1.err"
if [ -z $listfile ]; then
    listfile=/dev/stdout
    listerr=/dev/stderr
fi

source="/users/sheffler/data/decoys/"
dest='/users/sheffler/data/decoys/relax_native/'

#for i in `/bin/ls -F $source | grep /`; do
#    if [ -e $source/$i/gather.log ]; then
#	for j in `/bin/ls -F $source/$i | grep /`; do
#	    setup_relax_native $source $dest $i $j
#        done
#    fi
#done

i='phil_homolog'
for j in `/bin/ls -F $source/$i | grep /`; do
    setup_relax_native $source $dest $i $j
done


