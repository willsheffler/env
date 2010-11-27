

source="/users/sheffler/data/decoys/relax_native"
dest="/users/sheffler/data/decoys/"

tmpdir=`pwd`
for i in phil_homolog sraman_nmr kira_bench kira_hom bqian_caspbench; do
    #if [ -e $source/$i/gather.log ]; then
	for j in `/bin/ls -F $source/$i | grep /`; do
	    cd $dest/$i/$j
	    rm -rf $dest/$i/$j/relax_native*
	    #reorg_relax_native $source $dest $i $j
	    echo "cp $source/$i/$j/relax_native/ $dest/$i/$j/relax_native/"
	    mkdir -p $dest/$i/$j/relax_native/
	    cp $source/$i/$j/relax_native/* $dest/$i/$j/relax_native/
	    cp $source/$i/$j/rn${j:0:4}.fasc $dest/$i/$j/relax_native.fasc
	    ls relax_native/*.pdb > relax_native.list
        done
    #fi
done
cd $tmpdir
