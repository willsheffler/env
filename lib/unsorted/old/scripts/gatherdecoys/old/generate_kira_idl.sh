foreach i(1aa3 1ail 1bf4 1ig5 1r69 1tig 1vii 1acp 1b72 1erv 1pgx 1tif 1utg 256b)
    echo 'idealizing' $i
    cd ~/data/decoys/kira/$i
    ~/rosetta/rosetta_faidl.gcc -s rosetta_data/$i.pdb -idealize -fa_input > idealize.log
    cp ${i}_0001.pdb rosetta_data/
end
