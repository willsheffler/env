for set in bqian_caspbench kira_bench kira_hom phil_homolog sraman_nmr; do 
    for prot in `/bin/ls $set | grep -v gather`; do 
	#echo $set $prot
	wc -l $set/$prot/data/relax_native_region.data
	#for i in `/bin/ls $prot/data/*region.data`; do
	#    if [ ! -s $i ]; then
	#	echo removing $i
	#	rm -f $i
	#    fi
	#done
    done
done



