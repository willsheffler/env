source /users/sheffler/scripts/gatherdecoys/gatherdecoysfunc.sh

function gather {
    sns=/sraman/nmr/proteins/
    snt=/users/sheffler/data/decoys/sraman_nmr/
    loc1=(data      data       data        dump     dump     dump     dump    )
    loc2=(acyl/2abd aller/1bmw amyinh/2ait rnt/1ygw iif/1nha il8/1ikm il4/1bcn)
    prot=(1hb6      1who       1hoe        4rnt     1i27     3il8     2int    )
    for i in 0 1 2 3 4 5 6; do
	src=/${loc1[$i]}/$sns/${loc2[$i]}    
	gatherdecoys $src $snt/${prot[$i]} $src/native_trials/ng ${prot[$i]} \
	 	     _fa/ none _fa/ none $src/misc_files_1000
    done

    src=/dump/sraman/dsba/1a24/parsed_structures/parse_1/
    gatherdecoys $src $snt/1fvk $src/native_trials/ng 1fvk _fa/ none _fa/ none $src/misc_files_1000
}

snt=/users/sheffler/data/decoys/sraman_nmr/
mkdir -p $snt
gather > $snt/gather.log 2> $snt/gather.err


if [  ]; then

gatherdecoys ${sramandir}acyl/2abd/ac_fa/   ac1hb6.fasc ../misc_files_1000/ 1hb6
cp ${sramandir}acyl/2abd/native_trials/ng/ng*.pdb 1hb6/relax_native/

gatherdecoys ${sramandir}aller/1bmw/ac_fa/  ac1who.fasc ../misc_files_1000/ 1who
cp ${sramandir}aller/1bmw/native_trials/ng/ng*.pdb 1who/relax_native/

gatherdecoys ${sramandir}amyinh/2ait/ac_fa/ ac1hoe.fasc ../misc_files_1000/ 1hoe
cp ${sramandir}amyinh/2ait/native_trials/ng/ng*.pdb 1hoe/relax_native/

sramandir=/dump/sraman/nmr/proteins/

gatherdecoys ${sramandir}rnt/1ygw/ac_fa/    ac4rnt.fasc ../misc_files_1000/ 4rnt
cp ${sramandir}rnt/1ygw/native_trials/ng/ng*.pdb 4rnt/relax_native/

gatherdecoys ${sramandir}iif/1nha/ac_fa/    ac1i27.fasc ../misc_files_1000/ 1i27
cp ${sramandir}iif/1nha/native_trials/ng/acaa*.pdb 1i27/relax_native/

gatherdecoys ${sramandir}il8/1ikm/ad_fa/    ad3il8.fasc ../misc_files_1000/ 3il8
cp ${sramandir}il8/1ikm/native_trials/ng/ng*.pdb 3il8/relax_native/

gatherdecoys ${sramandir}il4/1bcn/ba_fa/    ba2int.fasc ../misc_files_1000/ 2int
cp ${sramandir}il4/1bcn/native_trials/ng/ng*.pdb 2int/relax_native/

gatherdecoys /dump/sraman/dsba/1a24/parsed_structures/parse_1/qa_fa/ qa1fvk.fasc ../misc_files_1000/ 1fvk
cp /dump/sraman/dsba/1a24/parsed_structures/parse_1/native_trials/ng/ng*.pdb 1fvk/relax_native/

fi
