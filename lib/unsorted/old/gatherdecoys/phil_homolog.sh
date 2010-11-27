
function copyphildecoys {
    source=$1
    dest=$2
    frags=$3
    prot=${4:0:4}
    chain=${4:4:1}
    suffix=${4:5:1}
    fragprefix=$5
    tmpdir=`pwd`
    echo gathering $prot $chain $suffix
    mkdir -p $dest/$prot
    cd $dest/$prot
    echo $source/$prot$chain$suffix > sourcedir.log
    cp -r $source/$prot$chain$suffix/rnd1/lowscore_decoys ./all
    /bin/ls all/*.pdb > all.list
    cp ~/data/decoys/paths.txt.template paths.txt
    python ~/scripts/setpdbchain.py '*' ' ' ./all/*.pdb > all_setchain.log
    cp -r $frags rosetta_data
    cd rosetta_data
    for j in `/bin/ls | egrep "^h001_"`; do
#	echo renaming $j
	cp $j ${prot}_${j:5}
    done
    if [ -e h001.pdb ]; then
	cp h001.pdb $prot.pdb
	~/scripts/setpdbchain.py '*' ' ' $prot.pdb > ${prot}_setchain.log
    elif [ -e h001_.pdb ]; then
	cp h001_.pdb $prot.pdb
	~/scripts/setpdbchain.py '*' ' ' $prot.pdb > ${prot}_setchain.log
    else
	echo "No pdb for" $prot $chain $suffix "!!!!"
    fi
    cp ${fragprefix}h001_03_05.200_v1_3 aa${prot}_03_05.200_v1_3
    cp ${fragprefix}h001_09_05.200_v1_3 aa${prot}_09_05.200_v1_3
    cd $tmpdir
}

function scoredecoys {
    location=$1
    prot=${2:0:4}
    chain=${2:4:1}
    tmpdir=`pwd`
    echo "scoring decoys in" $location/$prot
    cd $location/$prot
    nice ~/rosetta.gcc sc $prot _ -score -fa_input -l all.list -scorefile all > all_score.log 2> all_score.err
    cd $tmpdir
}

# source=/dump/pbradley/homolog_benchmark
# dest=/users/sheffler/data/decoys/phil_homolog
# prots="1af7_1 1di2a_ 1n0ua4 1dtja_ 1o2fb_ 
#        1mkya3 1ogwa_ 1dcja_ 1mla_2 2reb_2"

# for i in $prots; do    
#     copyphildecoys $source $dest /dump/pbradley/hom_fa/small_fold/frags/d$i/h001_ $i aa
# done

# copyphildecoys $source $dest /dump/pbradley/hom_fa/1b72_/h001_ 1b72A_ cc
# copyphildecoys $source $dest /dump/pbradley/hom_fa/1csp_/h001_ 1csp__ aa
# copyphildecoys $source $dest /dump/pbradley/hom_fa/1r69_/h001_ 1r69__ cc
# copyphildecoys $source $dest /dump/pbradley/hom_fa/1shfA/h001_ 1shfA_ aa
# copyphildecoys $source $dest /dump/pbradley/hom_fa/1tif_/h001_ 1tif__ aa
# copyphildecoys $source $dest /dump/pbradley/hom_fa/1tig_/h001_ 1tig__ aa

protswithnative="1af7_  1ogwa  1shfA  1tig_
                 1b72A  1dcja  1mla_  1r69_  1tif_"

#for i in $protswithnative; do
#    scoredecoys ~/data/decoys/phil_homolog $i
#done 

round2="1di2a_ 1n0ua4 1dtja_ 1o2fb_ 1mkya3 2reb_2"
for i in $round2; do
    scoredecoys ~/data/decoys/phil_homolog $i
done 

1dcj 1di2 1mla 1n0u 1ogw

prots="1af7 1b72 1csp 1dcj 1di2 1dtj 1mky 1mla 1n0u 1ogw 1r69 1shf 1tif 1tig 2reb"
for i in $prots; do
    echo $i    
#     mkdir -p $i/gsbr
#     mkdir -p $i/gsgr
#     head -n1 $i/all.fasc > $i/gsgr.fasc
#     head -n1 $i/all.fasc > $i/gsbr.fasc
#     sort -n -k32 $i/all.fasc | grep -v $i.pdb | grep -v filename | head -n 200 >> $i/gsgr.fasc
#     sort -n -k32 $i/all.fasc | grep -v $i.pdb | grep -v filename | tail -n 200 >> $i/gsbr.fasc
#     for j in `awk '{ print( $1) }' $i/gsgr.fasc | grep -v filename`; do
# 	cp $i/all/$j* $i/gsgr/
#     done
#     for j in `awk '{ print( $1) }' $i/gsbr.fasc | grep -v filename`; do
# 	cp $i/all/$j* $i/gsbr/
#     done
    cd $i
    /bin/ls gsgr/*.pdb > gsgr.list
    /bin/ls gsbr/*.pdb > gsbr.list
    cd ..
done

