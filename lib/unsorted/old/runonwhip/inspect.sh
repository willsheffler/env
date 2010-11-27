for set in bqian_caspbench kira_bench kira_hom phil_homolog sraman_nmr; do 
    for prot in `/bin/ls $set | grep -v gather`; do 
	for group in gsbr gsgr relax_native; do 
	    f=`ls $set/$prot/$group/*.decoyfeatures | wc -l`
	    p=`ls $set/$prot/$group/*.pdb | wc -l`
	    echo $set $prot $group $f $p
	#echo $set $prot
	#wc -l $set/$prot/data/relax_native_region.data
	#for i in `/bin/ls $prot/data/*region.data`; do
	#    if [ ! -s $i ]; then
	#	echo removing $i
	#	rm -f $i
	#    fi
	#done
	done
    done
done



