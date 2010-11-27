source /users/sheffler/scripts/gatherdecoys/gatherdecoysfunc.sh

bcs='/dump/bqian/data/bqian/projects/caspbench/'
bct='/users/sheffler/data/decoys/bqian_caspbench/'
bcnat='/users/sheffler/data/relax_native/bqian_caspbench'
bcp='196 199 205 234'
bcdt='decoys_r'
bcdl='none'
bcst='decoys_r'
bcsl='none'

function gather {
for i in $bcp; do
 gatherdecoys $bcs/CMHARD/T0$i $bct/t$i $bcnat/t$i/relax_native t$i $bcdt $bcdl $bcst $bcsl $bcs/CMHARD/T0$i
done

bcp='206 223 224'
for i in $bcp; do
 gatherdecoys $bcs/FRHOM/T0$i  $bct/t$i $bcnat/t$i/relax_native t$i $bcdt $bcdl $bcst $bcsl $bcs/FRHOM/T0$i
done

gatherdecoys $bcs/CMHARD/T2321 $bct/t232 $bcnat/t232/relax_native t232 $bcdt $bcdl $bcst $bcsl $bcs/CMHARD/T2321
gatherdecoys $bcs/FRHOM/T02492 $bct/t249 $bcnat/t249/relax_native t249 $bcdt $bcdl $bcst $bcsl $bcs/FRHOM/T02492
}

mkdir -p '/users/sheffler/data/decoys/bqian_caspbench/'
gather > $bct/gather.log 2> $bct/gather.err
