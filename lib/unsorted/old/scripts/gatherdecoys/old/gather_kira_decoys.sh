#!/bin/tcsh

foreach i (1aa3 1acp 1ail 1b72 1bf4 1erv 1ig5 1gx 1r69 1tif 1tig 1utg 1vii 256b)
    echo 'gathering kira decoys' $i
    cd ~/data/decoys/kira/
    #mkdir $i
    cd $i

    cp -r /data/kira/fullatom_benchmark/$i/lowscore lowscore_decoys
    foreach j (`du ./lowscore_decoys/* | grep '^0' | awk '{print $2}'`)
	echo 'removing' $i $j
	rm -f $j
    end
    /bin/ls lowscore_decoys/*.pdb > lowscore.list

    cp -r /data/kira/fullatom_benchmark/$i/lowrms lowrms_decoys
    foreach j (`du ./lowrms_decoys/* | grep '^0' | awk '{print $2}'`)
	echo 'removing' $i $j
	rm -f $j
    end
    /bin/ls lowrms_decoys/*.pdb > lowrms.list

    mkdir gsbr_decoys
    touch gsbr.list


    #/bin/ls native/*.pdb > native.list
    #mkdir rosetta_data
    #cp /data/kira/fullatom_benchmark/fragments/*$i* rosetta_data
    #cp /data/sheffler/decoys/paths.txt.template paths.txt
end
