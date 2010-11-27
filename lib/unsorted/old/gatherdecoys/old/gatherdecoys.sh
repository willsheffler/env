#~/bin/bash

function gatherdecoys 
{
    decoydir=$1
    scorefile=$2
    rosettadir=$3
    prot=$4

    echo gathering decoys from $decoydir
    echo using scorefile $scorefile
    echo putting decoys into dir $prot
    mkdir ${prot}
    cd ${prot}
    mkdir lowscore_decoys
    mkdir lowrms_decoys
    mkdir gsbr_decoys
    mkdir relax_native
    mkdir rosetta_data
    mkdir tmp

    cp ${decoydir}$scorefile all.fasc

    head -n1 all.fasc > lowscore.fasc
    sort -k2 -n all.fasc | grep -av ${prot}.pdb | grep -av filename | head -n200 >> lowscore.fasc
    head -n1 all.fasc > lowrms.fasc
    sort -k32 -n all.fasc | grep -av ${prot}.pdb | grep -av filename | head -n200 >> lowrms.fasc

    ndecoy=`wc -l all.fasc | awk '{print $1}'`
    tmp=`expr $ndecoy / 5`
    nrms=`expr $tmp \* 4`
    sort -k32 -n all.fasc | grep -av ${prot}.pdb | grep -av filename | tail -n$nrms > tmp/rmsq20.fasc
    head -n1 all.fasc > gsbr.fasc
    sort -k2 -n tmp/rmsq20.fasc | grep -av filename | head -n200 >> gsbr.fasc
    #rm -f tmp/rmsq20.fasc

    nice R CMD BATCH ~/scripts/gatherdecoys/gather_decoys_plot_scores.R
    #rm gather_decoys_plot_scores.Rout

    ~/scripts/gatherdecoys/copypdbs.py lowscore.fasc $decoydir lowscore_decoys
    ~/scripts/gatherdecoys/copypdbs.py lowrms.fasc $decoydir lowrms_decoys
    ~/scripts/gatherdecoys/copypdbs.py gsbr.fasc $decoydir gsbr_decoys

    /bin/ls lowscore_decoys/*.pdb > lowscore.list
    /bin/ls lowrms_decoys/*.pdb > lowrms.list
    /bin/ls gsbr_decoys/*.pdb > gsbr.list

    cp ${decoydir}${rosettadir}*${prot}* rosetta_data/
    
    cd ..
}

decoydir=${bindir}CMHARD/T0196/decoys_r2/
scorefile=r2t196.fasc
rosettadir=../
prot=t196


destdir=~/data/decoys/bqian/
cd $destdir
bindir=/users/bqian/data/bqian/projects/caspbench/

gatherdecoys ${bindir}CMHARD/T0196/decoys_r2/ r2t196.fasc ../ t196
gatherdecoys ${bindir}CMHARD/T0199/decoys_r2/ r2t199.fasc ../ t199
gatherdecoys ${bindir}CMHARD/T0205/decoys_r2/ r2t205.fasc ../ t205
gatherdecoys ${bindir}CMHARD/T0223/decoys_r1/ r1t223.fasc ../ t223
gatherdecoys ${bindir}CMHARD/T2321/decoys_r2/ r2t232.fasc ../ t232
gatherdecoys ${bindir}CMHARD/T0234/decoys_r2/ r2t234.fasc ../ t234
gatherdecoys ${bindir}CMHARD/T0279/decoys_r1/ r1t279.fasc ../ t279
gatherdecoys ${bindir}FRHOM/T0206/decoys_r5/ r5t206.fasc ../ t206
gatherdecoys ${bindir}FRHOM/T0249/decoys_r5/ r5t249.fasc ../ t249
gatherdecoys ${bindir}FRHOM/T0262/decoys_r5/ r5t262.fasc ../ t262
rm t262/rosetta_data/rbt262_dump00*



destdir=~/data/decoys/sraman/
cd $destdir
sramandir=/data/sraman/nmr/proteins/

gatherdecoys ${sramandir}acyl/2abd/ac_fa/   ac1hb6.fasc ../misc_files_1000/ 1hb6
cp ${sramandir}acyl/2abd/native_trials/ng/ng*.pdb 1hb6/relax_native/
cp ${sramandir}acyl/2abd/native_trials/ng/ng1hb6.fasc 1hb6/relax_native.fasc
cd 1hb6; /bin/ls relax_native/*.pdb > relax_native.list; cd ..

gatherdecoys ${sramandir}aller/1bmw/ac_fa/  ac1who.fasc ../misc_files_1000/ 1who
cp ${sramandir}aller/1bmw/native_trials/ng/ng*.pdb 1who/relax_native/
cp ${sramandir}aller/1bmw/native_trials/ng/ng1who.fasc 1who/relax_native.fasc
cd 1who; /bin/ls relax_native/*.pdb > relax_native.list; cd ..

gatherdecoys ${sramandir}amyinh/2ait/ac_fa/ ac1hoe.fasc ../misc_files_1000/ 1hoe
cp ${sramandir}amyinh/2ait/native_trials/ng/ng*.pdb 1hoe/relax_native/
cp ${sramandir}amyinh/2ait/native_trials/ng/ng1hoe.fasc 1hoe/relax_native.fasc
cd 1hoe; /bin/ls relax_native/*.pdb > relax_native.list; cd ..

sramandir=/dump/sraman/nmr/proteins/

gatherdecoys ${sramandir}rnt/1ygw/ac_fa/    ac4rnt.fasc ../misc_files_1000/ 4rnt
cp ${sramandir}rnt/1ygw/native_trials/ng/ng*.pdb 4rnt/relax_native/
cp ${sramandir}rnt/1ygw/native_trials/ng/ng4rnt.fasc 4rnt/relax_native.fasc
cd 4rnt; /bin/ls relax_native/*.pdb > relax_native.list; cd ..

gatherdecoys ${sramandir}iif/1nha/ac_fa/    ac1i27.fasc ../misc_files_1000/ 1i27
cp ${sramandir}iif/1nha/native_trials/ng/acaa*.pdb 1i27/relax_native/
cp ${sramandir}iif/1nha/native_trials/ng/ac1i27.fasc 1i27/relax_native.fasc
cd 1i27; /bin/ls relax_native/*.pdb > relax_native.list; cd ..

gatherdecoys ${sramandir}il8/1ikm/ad_fa/    ad3il8.fasc ../misc_files_1000/ 3il8
cp ${sramandir}il8/1ikm/native_trials/ng/ng*.pdb 3il8/relax_native/
cp ${sramandir}il8/1ikm/native_trials/ng/ng3il8.fasc 3il8/relax_native.fasc
cd 3il8; /bin/ls relax_native/*.pdb > relax_native.list; cd ..

gatherdecoys ${sramandir}il4/1bcn/ba_fa/    ba2int.fasc ../misc_files_1000/ 2int
#cp ${sramandir}il4/1bcn/native_trials/ng/ng*.pdb 2int/relax_native/
#cp ${sramandir}il4/1bcn/native_trials/ng/ng2int.fasc 2int/relax_native.fasc
#cd 2int; /bin/ls relax_native/*.pdb > relax_native.list; cd ..

gatherdecoys /dump/sraman/dsba/1a24/parsed_structures/parse_1/qa_fa/ qa1fvk.fasc ../misc_files_1000/ 1fvk
cp /dump/sraman/dsba/1a24/parsed_structures/parse_1/native_trials/ng/ng*.pdb 1fvk/relax_native/
cp /dump/sraman/dsba/1a24/parsed_structures/parse_1/native_trials/ng/ng1fvk.fasc 1fvk/relax_native.fasc
cd 1fvk; /bin/ls relax_native/*.pdb > relax_native.list; cd ..
