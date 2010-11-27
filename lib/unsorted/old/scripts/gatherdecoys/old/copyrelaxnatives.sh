function cprn {
    name=$1
    prot=$2
    echo $name $prot
    src=/users/sheffler/data/decoys/$name/$prot/
    tgt=/users/sheffler/data/relax_native/$name/$prot
    mkdir -p $tgt/relax_native/$name/$i
    cp -r $src/relax_native* $tgt
}


name=bqian
for i in t196  t199  t205  t206  t223  t232  t234  t249  t262  t279; do
    cprn $name $i
done

name=kira
for i in 1aa3 1acp 1ail 1b72 1bf4 1erv 1ig5 1pgx 1r69 1tif 1tig 1utg 1vii 256b; do
    cprn $name $i
done

name=sraman
for i in 1fvk  1hb6  1hoe  1who  3il8  4rnt; do
    cprn $name $i
done
