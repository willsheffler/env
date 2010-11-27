#!/bin/bash

dir=~/data/decoys/bqian/t196/
list=tmp.list
log=~/decoy_buried_polar/test.bp.log
prot=t196

whips=(whip01 whip02 whip03 whip04 whip05 whip06 whip07 whip08 whip09 whip10 whip11)
whip=0

function rundecoystats
{
    dir=$1; list=$2; log=$3; prot=$4

    args='-decoyfeatures -score -fa_input -try_both_his_tautomers -constant_seed'
    cmd="nice ~/rosetta_features.gcc bp $prot _ ${args} -l ${list}"
    ssh ${whips[whip]} "cd $dir; rm -f bp${prot}.fasc; $cmd > ${log}; " &
    tmp=`expr $whip + 1`
    whip=`expr $tmp % 11`
}

#rundecoystats ~/data/decoys/bqian/t196/ tmp.list ~/tmp.log t196

dirh=~/data/decoys/
logh=~/decoy_buried_polar/log/

name=bqian
for i in t196  t199  t205  t206  t223  t232  t234  t249  t262  t279; do
    dir="~/data/decoys/${name}/$i"
    rundecoystats $dir gsbr.list         ${logh}gsbr.$name.$i.bp.log $i
    rundecoystats $dir relax_native.list ${logh}relax_native.$name.$i.bp.log $i
done

name=kira
for i in 1aa3 1acp 1ail 1b72 1bf4 1erv 1ig5 1pgx 1r69 1tif 1tig 1utg 1vii 256b; do
    dir="~/data/decoys/${name}/$i"
    rundecoystats $dir lowscore.list     ${logh}lowscore.$name.$i.bp.log $i
    rundecoystats $dir relax_native.list ${logh}relax_native.$name.$i.bp.log $i
done

name=sraman
for i in 1fvk  1hb6  1hoe  1who  3il8  4rnt; do
    dir="~/data/decoys/${name}/$i"
    rundecoystats $dir gsbr.list     ${logh}gsbr.$name.$i.bp.log $i
    rundecoystats $dir relax_native.list ${logh}relax_native.$name.$i.bp.log $i
done



for i in 01 02 03 04 05 06 07 08 09 10 11
do
  echo $i 
  ssh whip$i killall rosetta_features.gcc
  #ssh whip$i ps aux | grep sheffler | grep rosetta
done

tot=0
for i in 01 02 03 04 05 06 07 08 09 10 11
do
  tmp=`ssh whip$i 'ps aux | grep rosetta_features | grep -v grep | wc -l'`
  echo $i $tmp
  tot=`expr $tot + $tmp`
done
echo $tot

cd ~/decoy_buried_polar/
for i in `/bin/ls log`; do
    grep DF_POLAR log/$i > df_polar/$i.dfpolar    
done

rm -f relax_native_buried_polar.data decoy_bruied_polar.data

for i in `/bin/ls df_polar/*native*`; do
    cat $i >> relax_native_buried_polar.data
done

for i in `/bin/ls df_polar/*lowscore*`; do
    cat $i >> decoy_bruied_polar.data
done

for i in `/bin/ls df_polar/*gsbr*`; do
    cat $i >> decoy_bruied_polar.data
done

###################3
# pdbs
#  -decoyfeatures -score -fa_input -try_both_his_tautomers -chain_inc -constant_seed -skip_missing_residues -l pdb.list6
###################

# do on whip
cd ~/pdb_buried_polar
split -d -a1 -l 400 ~/data/pdb/cullpdb_pc50_res2.0_R0.25_d050604_chains3866.list pdb.list


whips=(whip01 whip02 whip03 whip04 whip05 whip06 whip07 whip08 whip09 whip10 whip11)
whip=0
function runpdb
{
    dir=$1; list=$2; log=$3;
    args='-decoyfeatures -score -fa_input -try_both_his_tautomers -chain_inc -constant_seed'
    cmd="nice ~/rosetta_features.gcc ${args} -l ${list}"
    ssh ${whips[whip]} "cd $dir; $cmd > ${log}; " &
    tmp=`expr $whip + 1`
    whip=`expr $tmp % 11`
}

for i in 0 1 2 3 4 5 6 7 8 9; do
    runpdb ~/pdb_buried_polar/ pdb.list$i pdb_buried_polar_$i.log
done

for i in 0 1 2 3 4 5 6 7 8 9; do
    cat pdb_buried_polar_$i.log >> pdb_buried_polar.log
done

