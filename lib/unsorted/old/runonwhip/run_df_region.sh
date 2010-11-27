
###################3
# pdbs
#  -decoyfeatures -score -fa_input -try_both_his_tautomers -chain_inc -constant_seed -skip_missing_residues -l pdb.list6
###################


whips=(whip01 whip02 whip03 whip04 whip05 whip06 whip08 whip09 whip10 whip11)
whip=1
function runpdb
{
    dir=$1; list=$2; log=$3; job=$4
    args="-decoyfeatures $job -score -fa_input -try_both_his_tautomers  -constant_seed"
    cmd="nice -n15 ~/rosetta/rosetta_features/rosetta.gcc ${args} -l ${list}"
    ssh -x ${whips[whip]} "cd $dir; $cmd > $log.log 2> $log.err; " &
    tmp=`expr $whip + 1`
    whip=`expr $tmp % ${#whips[*]}`
}

whip=0
for i in 0 1 2 3 4 5 6 7 8 9; do
    runpdb ~/pdb_energies/ pdb.list$i pdb_region_$i region
done

#for i in 0 1 2 3 4 5 6 7 8 9; do
#    cat pdb_energies_$i.log >> pdb_buried_polar.log
#done

#for i in `/bin/ls | grep -v gather`; do echo $i `grep 'output features score' $i/log/all_score.log | wc -l` `grep 'output features score' $i/log/gsbr_score.log | wc -l` `grep 'output features region' $i/log/gsgr_region.log | wc -l` 'output features region' $i/log/all_region.log | wc -l` `grep 'output features region' $i/log/gsbr_region.log | wc -l` `grep 'output features region' $i/log/gsgr_region.log | wc -l` 'output features region' $i/log/all_region.log | wc -l` `grep 'output features region' $i/log/gsbr_region.log | wc -l` `grep 'output features region' $i/log/gsgr_region.log | wc -l`; done
