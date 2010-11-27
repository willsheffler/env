function removenullpdbs {
    dir=$1
    for i in `/bin/ls $dir/*.pdb `; do
	#echo $i
	if [ ! -s $i ]; then
	    echo rm -f $i
	    rm -f $i
	fi
    done
}
