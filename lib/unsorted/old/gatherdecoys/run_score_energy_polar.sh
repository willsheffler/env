function run_score_energy_polar {

    dir=$1
    set=$2
    prot=$3
    tmpdir=`pwd`
    echo $dir $set $prot
    cd $dir/$set/$prot
    chain=`/bin/ls rosetta_data | egrep "..$j.03......._v..3"`
    chain=${chain:6:1}
    cp ~/data/decoys/paths.txt.template paths.txt
    mkdir -p log
    rm -f *_score.fasc
    rm -f *.log
    rm -f *.data
    rm -f *.err
    for job in gsgr gsbr all relax_native; do # relax_native `/bin/ls -d rand_*` `/bin/ls -d rms*` lowscore lowrms; do
	if [ -d $job ]; then
	    echo $dir $set $prot $job
#	    if [ ! -s $job.list ]; then
#		echo removing null dir $job
#		rmdir $job
#		rm -f $job.*
#		continue
#	    fi
#	    if [ -e ${job}_energy.data ]; then
#		echo skipping $job, already done
#		continue
#	    fi
#	    for i in `/bin/ls $job/*.pdb`; do
#		if [ ! -e ${i:0:`expr length $i - 4`}.dssp ]; then
#		    nice -n15 ~pbradley/dssp $i > ${i:0:`expr length $i - 4`}.dssp
#		fi
#	    done
 	    #rm -f $job/$prot.pdb
 	    #/bin/ls $job/*.pdb > $job.list
	    #
 	    #nice -n15 ~/rosetta/rosetta_features/rosetta.gcc sc $prot $chain -score -scorefile ${job}_score -fa_input -try_both_his_tautomers -l $job.list > log/${job}_score.log 2> log/${job}_score.err
 	    #nice -n15 ~/rosetta/rosetta_features/rosetta.gcc -score ${job} -decoyfeatures energyandpolar -fa_input -try_both_his_tautomers -l $job.list > log/${job}_energy_polar.log 2> log/${job}_energy_polar.err	   
 	    #nice -n15 grep DF_POLAR log/${job}_energy_polar.log > ${job}_polar.data
 	    #nice -n15 grep DF_ENERG log/${job}_energy_polar.log > ${job}_energy.data
 	    #if [ ! -s ${job}_polar.data ] ; then
 	#	rm -f ${job}_polar.data ${job}_energy.data
 	#    fi
 	#    if [ ! -s ${job}_energy.data ] ; then
 	#	rm -f ${job}_polar.data ${job}_energy.data
	    #    fi
	    rm -f log/tmp_score_file.fasc
	    rm data/${job}_region.data
	    nice -n15 ~/rosetta/rosetta_features/rosetta.gcc rg $prot $chain -score ${job} -decoyfeatures region -fa_input -try_both_his_tautomers -l $job.list > log/${job}_region.log 2> log/${job}_region.err	-scorefile log/tmp_score_file   
	    nice -n15 grep DF_REGION log/${job}_region.log > data/${job}_region.data
	fi
    done
    cd $tmpdir

}

run_score_energy_polar $1 $2 $3
