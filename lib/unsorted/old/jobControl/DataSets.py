prots = {}
prots['kira_bench']      = ['1acp', 
                            '1b72',   
                            '1ig5', 
                            '1r69', 
                            '1tig', 
                            '1vii', 
                            '1aa3', 
                            '1ail', 
                            '1bf4', 
                            '1pgx', 
                            '1tif', 
                            '1utg', 
                            '256b']
prots['phil_homolog']    = ['1af7_', 
                            '1csp_', 
                            '1di2_', 
                            '1mky_', 
                            '1n0u_', 
                            '1ogw_', 
                            '1shf_', 
                            '1tig_', 
                            '1b72_', 
                            '1dcj_', 
                            '1dtj_', 
                            '1mla_', 
                            '1o2f_', 
                            '1r69_', 
                            '1tif_', 
                            '2reb_']
prots['kira_hom']        = ['1a4pA', 
                            '1agi', 
                            '1ar0', 
                            '1b07', 
                            '1bmb', 
                            '1cfy', 
                            '1d2n', 
                            '1erv', 
                            '1pva', 
                            '1acf_', 
                            '1ahn', 
                            '1aw2', 
                            '1be9', 
                            '1btn', 
                            '1crb', 
                            '1dol', 
                            '1ig5', 
                            '1qav']
prots['bqian_caspbench'] = ['t196_', 
                            't199_', 
                            't205_', 
                            't206_', 
                            't223_', 
                            't224_', 
                            't232_', 
                            't234_', 
                            't249_']
prots['sraman_nmr']      = ['1fvk', 
                            '1hb6', 
                            '1hoe', 
                            '1i27', 
                            '1who', 
                            '2int', 
                            '3il8', 
                            '4rnt']

for k in prots:
    prots[k].sort()

for k in prots:
    print 'dataSet["'+k+'"] = '+" "*(20-len(k))+'[',
    for p in prots[k]:
        if len(p) is 4:
            print '"'+p+'_",'
        else: print '"'+p+'",'
    print ']'
