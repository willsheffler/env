foreach i (`/bin/ls`)
    echo `ls $i/relax_native/*.pdb | wc -l`  `du -s $i/relax_native/*.pdb | sort -n -k1 | head -n1`
end
