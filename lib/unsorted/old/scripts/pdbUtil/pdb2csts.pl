#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.8 $
##  $Date: 2004/11/18 00:07:50 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering

$def_cst_radius = 10.0;  # in angstroms
$def_cst_pad    =  2.0;  # in angstroms
$def_cst_isep   =    5;  # in res, so i can compare with j >= i+5

#$cst_radius_start       =  4.0;  # in angstroms 
$cst_radius_stop         = 20.0;  # in angstroms 
$cst_radius_incr         =  2.0;  # in angstroms
$long_range_isep         =   15;  # in res
#$long_range_cst_num_stop =    6;  # min num of long range constraints
$long_range_cst_num_stop =    3;  # min num of long range constraints

###############################################################################
# init
###############################################################################

my %opts = &getCommandLineOptions ();
$pdb                 = $opts{pdb};
$cst_file            = $opts{cstout};
$cst_radius          = $opts{radius};
$cst_pad             = $opts{pad};
$cst_isep            = $opts{isep};
$cstCA_flag          = $opts{cstCA};
$cstCB_flag          = $opts{cstCB};
$add_long_range_flag = $opts{addlongrange};
$add_long_range_flag = $opts{addlongrange};

###############################################################################
# main
###############################################################################

# read pdb for occupancy
#
@pdbbuf = &fileBufArray ($pdb);
$last_res_n = undef;
@fasta = ();
@occ   = ();
for ($i=0; $i <= $#pdbbuf; ++$i) {
    $line = $pdbbuf[$i];
    if ($line =~ /^ATOM/) {

	$atomname = substr ($line, 12, 4);
	if ($atomname eq ' CA ') {
	    $res3 = substr ($line, 17, 3);
	    push (@fasta, &mapResCode($res3));
	}

	$atom_occ = substr ($line, 54, 6);
	next if ($atom_occ <= 0.0);

	$res_n = substr ($line, 22, 4);
	$res_i = $res_n - 1;
	$last_res_n = $res_n  if (! defined $last_res_n);
	if ($res_n != $last_res_n || $i == $#pdbbuf) {
	    if ($i == $#pdbbuf) {
		$atomname = substr ($line, 12, 4);
		if ($atomname eq ' N  ') {
		    $N_occ  = 1;
		} elsif ($atomname eq ' CA ') {
		    $CA_occ = 1;
		} elsif ($atomname eq ' C  ') {
		    $C_occ  = 1;
		} elsif ($atomname eq ' O  ') {
		    $O_occ  = 1;
		}
		$occ[$res_i] = 1  if ($N_occ && $CA_occ && $C_occ && $O_occ);
	    }
	    $last_res_n = $res_n;
	    $N_occ      = undef;
	    $CA_occ     = undef;
	    $C_occ      = undef;
	    $O_occ      = undef;
	}
	$atomname = substr ($line, 12, 4);
	if ($atomname eq ' N  ') {
	    $N_occ  = 1;
	} elsif ($atomname eq ' CA ') {
	    $CA_occ = 1;
	} elsif ($atomname eq ' C  ') {
	    $C_occ  = 1;
	} elsif ($atomname eq ' O  ') {
	    $O_occ  = 1;
	}
	$occ[$res_i] = 1  if ($N_occ && $CA_occ && $C_occ && $O_occ);
    }
}
$seq_len = $#fasta + 1;
# debug
#for ($i=0; $i < $seq_len; ++$i) {
#    print "occ[$i] = '$occ[$i]'\n";
#}


# read pdb from which to measure squared distances
#
foreach $line (@pdbbuf) {
    if ($line =~ /^ATOM/) {
	$atomtype          = substr ($line, 12, 4);
	$restype           = substr ($line, 17, 3);
	$res_i             = substr ($line, 22, 4) - 1;

	# N, CA, C, O, CB only (CB recalculated below)
	if ($occ[$res_i] > 0 && defined &typeI($atomtype)) {
	    $src->[$res_i]->[&typeI($atomtype)]->[0] = substr ($line, 30, 8);
	    $src->[$res_i]->[&typeI($atomtype)]->[1] = substr ($line, 38, 8);
	    $src->[$res_i]->[&typeI($atomtype)]->[2] = substr ($line, 46, 8);
	}
	if (!  defined $src->[$res_i]->[&typeI('CB')]
	    && defined $src->[$res_i]->[&typeI('N')]
	    && defined $src->[$res_i]->[&typeI('CA')]
	    && defined $src->[$res_i]->[&typeI('C')]
	    && $occ[$res_i] > 0
	    ) {
	    $src->[$res_i]->[&typeI('CB')] = 
		&getCbCoords ($src->[$res_i]->[&typeI('N')],
			      $src->[$res_i]->[&typeI('CA')],
			      $src->[$res_i]->[&typeI('C')],
			      $restype);
	}
    }
}


# measure pdb distances (squared for efficiency)
#
if ($cstCA_flag) {
    for ($q_i=0; $q_i <= $#{$src}; ++$q_i) {
	$p_dist_sq->[$q_i]->[$q_i]->{CA}->{CA} = 0;
	for ($q_j=$q_i+1; $q_j <= $#{$src}; ++$q_j) {
	    if ($occ[$q_i] > 0 && $occ[$q_j] > 0) {
		$p_dist_sq->[$q_i]->[$q_j]->{CA}->{CA}
		= $p_dist_sq->[$q_j]->[$q_i]->{CA}->{CA} = &measureDistSq ($src->[$q_i], 
									   $src->[$q_j],
									   'CA',
									   'CA');
	    }
	}
    }
}
if ($cstCB_flag) {
    for ($q_i=0; $q_i <= $#{$src}; ++$q_i) {
	$p_dist_sq->[$q_i]->[$q_i]->{CB}->{CB} = 0;
	for ($q_j=$q_i+1; $q_j <= $#{$src}; ++$q_j) {
	    if ($occ[$q_i] > 0 && $occ[$q_j] > 0) {
		$p_dist_sq->[$q_i]->[$q_j]->{CB}->{CB}
		= $p_dist_sq->[$q_j]->[$q_i]->{CB}->{CB} = &measureDistSq ($src->[$q_i], 
									   $src->[$q_j],
									   'CB',
									   'CB');
	    }
	}
    }
}


# gradually loosen radius
#   stop when enough long range constraints in multiple elements
#
$CA_csts = +[];
$CB_csts = +[];
@long_range_constraint_cnt = ();
$long_range_csts = +[];

$stop_radius = ($add_long_range_flag) ? $cst_radius_stop : $cst_radius;

for ($cst_radius_iter  = $cst_radius; 
     $cst_radius_iter <= $stop_radius; 
     $cst_radius_iter += $cst_radius_incr
     ) {
    
#    print "CST_RADIUS_ITER: '$cst_radius_iter'\n";  # DEBUG
    
    $cst_radius_iter_sq = $cst_radius_iter * $cst_radius_iter;
    
    @long_range_constraint_cnt_this_round = ();
    
    for ($q_i=0; $q_i <= $#fasta; ++$q_i) {
	for ($q_j=$q_i+$cst_isep; $q_j <= $#fasta; ++$q_j) {
	    
	    if ($cst_radius_iter == $cst_radius ||  # accept all first pass
		($q_j - $q_i >= $long_range_isep &&
		 ($long_range_constraint_cnt[$q_i] < $long_range_cst_num_stop ||
		  $long_range_constraint_cnt[$q_j] < $long_range_cst_num_stop)
		 )
		) {
		
		if ($occ[$q_i] > 0 && $occ[$q_j] > 0) {
		    if ($cstCA_flag) {
			if (! defined $CA_csts->[$q_i]->[$q_j]) {
			    $dist_sq = $p_dist_sq->[$q_i]->[$q_j]->{CA}->{CA};
			    if ($dist_sq <= $cst_radius_iter_sq) {
				$dist = sqrt ($dist_sq);
				
				if ($cst_radius_iter == $cst_radius) {
				    $CA_csts->[$q_i]->[$q_j] = $dist;
				}
				if ($q_j - $q_i >= $long_range_isep) {
				    $long_range_csts->[$q_i] = +[]  if (! defined $long_range_csts->[$q_i]);
				    $long_range_csts->[$q_j] = +[]  if (! defined $long_range_csts->[$q_j]);
				    
				    # only count diverse long range constraints (themselves well separated)
				    $new_long_range_cst = 'true';
				    for ($lr_i=0; $lr_i <= $#{$long_range_csts->[$q_i]}; ++$lr_i) {
					$q_j2 = $long_range_csts->[$q_i]->[$lr_i];
					next if ($q_j2 == $q_j);
					if (abs ($q_j2 - $q_j) < $long_range_isep) {
					    $new_long_range_cst = undef;
					    last;
					}
				    }
				    if ($new_long_range_cst) {
					++$long_range_constraint_cnt_this_round[$q_i];
					push (@{$long_range_csts->[$q_i]}, $q_j);
					push (@{$long_range_csts->[$q_j]}, $q_i);
					$CA_csts->[$q_i]->[$q_j] = $dist;
				    }
				    $new_long_range_cst = 'true';
				    for ($lr_i=0; $lr_i <= $#{$long_range_csts->[$q_j]}; ++$lr_i) {
					$q_i2 = $long_range_csts->[$q_j]->[$lr_i];
					next if ($q_i2 == $q_i);
					if (abs ($q_i2 - $q_i) < $long_range_isep) {
					    $new_long_range_cst = undef;
					    last;
					}
				    }
				    if ($new_long_range_cst) {
					++$long_range_constraint_cnt_this_round[$q_j];
					push (@{$long_range_csts->[$q_i]}, $q_j);
					push (@{$long_range_csts->[$q_j]}, $q_i);
					$CA_csts->[$q_i]->[$q_j] = $dist;
				    }
				}					    
			    }
			}
		    }
		    if ($cstCB_flag && 
			$fasta[$q_i] !~ /G/i &&           # no glycine CB
			$fasta[$q_j] !~ /G/i
			) {
			if (! defined $CB_csts->[$q_i]->[$q_j]) {
			    $dist_sq = $p_dist_sq->[$q_i]->[$q_j]->{CB}->{CB};
			    if ($dist_sq <= $cst_radius_iter_sq) {
				$dist = sqrt ($dist_sq);
				
				if ($cst_radius_iter == $cst_radius) {
				    $CB_csts->[$q_i]->[$q_j] = $dist;
				}
				if ($q_j - $q_i >= $long_range_isep) {
				    $long_range_csts->[$q_i] = +[]  if (! defined $long_range_csts->[$q_i]);
				    $long_range_csts->[$q_j] = +[]  if (! defined $long_range_csts->[$q_j]);
				    
				    # only count diverse long range constraints (themselves well separated)
				    $new_long_range_cst = 'true';
				    for ($lr_i=0; $lr_i <= $#{$long_range_csts->[$q_i]}; ++$lr_i) {
					$q_j2 = $long_range_csts->[$q_i]->[$lr_i];
					next if ($q_j2 == $q_j);
					if (abs ($q_j2 - $q_j) < $long_range_isep) {
					    $new_long_range_cst = undef;
					    last;
					}
				    }
				    if ($new_long_range_cst) {
					++$long_range_constraint_cnt_this_round[$q_i];
					push (@{$long_range_csts->[$q_i]}, $q_j);
					push (@{$long_range_csts->[$q_j]}, $q_i);
					$CB_csts->[$q_i]->[$q_j] = $dist;
				    }
				    $new_long_range_cst = 'true';
				    for ($lr_i=0; $lr_i <= $#{$long_range_csts->[$q_j]}; ++$lr_i) {
					$q_i2 = $long_range_csts->[$q_j]->[$lr_i];
					next if ($q_i2 == $q_i);
					if (abs ($q_i2 - $q_i) < $long_range_isep) {
					    $new_long_range_cst = undef;
					    last;
					}
				    }
				    if ($new_long_range_cst) {
					++$long_range_constraint_cnt_this_round[$q_j];
					push (@{$long_range_csts->[$q_i]}, $q_j);
					push (@{$long_range_csts->[$q_j]}, $q_i);
					$CB_csts->[$q_i]->[$q_j] = $dist;
				    }
				}					    
			    }
			}
		    }
		}
	    }
	}
    }
    
    # add in new constraint counts
    #
    for ($q_i=0; $q_i <= $#fasta; ++$q_i) {
	$long_range_constraint_cnt[$q_i] += $long_range_constraint_cnt_this_round[$q_i];
#    print "ADDED NEW LONG RANGE CNTS for res $q_i: '".$long_range_constraint_cnt_this_round[$q_i] ."' -> '". $long_range_constraint_cnt[$q_i] ."'\n";  # DEBUG
    }
}
    

# put into constraints format
#
@constraints = ();
for ($q_i=0; $q_i <= $#fasta; ++$q_i) {
    if (defined $CA_csts->[$q_i]) {
	for ($q_j=$q_i+$cst_isep; $q_j <= $#fasta; ++$q_j) {
	    if (defined $CA_csts->[$q_i]->[$q_j]) {
		
		$dist = $CA_csts->[$q_i]->[$q_j];
		$max_dist = $dist + $cst_pad;
		$min_dist = $dist - $cst_pad;
		$true_dist = 0.0;
		push (@constraints, 
		      sprintf (" %6d %2s %6d %2s   %10.2f %10.2f %10.2f",
			       $q_i+1, 'CA', $q_j+1, 'CA', $max_dist, $min_dist, $true_dist));
	    }
	}
    }
    if (defined $CB_csts->[$q_i]) {
	for ($q_j=$q_i+$cst_isep; $q_j <= $#fasta; ++$q_j) {
	    if (defined $CB_csts->[$q_i]->[$q_j]) {
		
		$dist = $CB_csts->[$q_i]->[$q_j];
		$max_dist = $dist + $cst_pad;
		$min_dist = $dist - $cst_pad;
		$true_dist = 0.0;
		push (@constraints, 
		      sprintf (" %6d %2s %6d %2s   %10.2f %10.2f %10.2f",
			       $q_i+1, 'CB', $q_j+1, 'CB', $max_dist, $min_dist, $true_dist));
	    }
	}
    }
}


# write cst file
#
$num_recs = $#constraints + 1;
@header = ("NMR_v3.0");
if ($cstCA_flag && $cstCB_flag) {
    push (@header, "CA-CA and CB-CB csts from ksync alignment");
} elsif ($cstCA_flag) {
    push (@header, "CA-CA csts from ksync alignment");
} elsif ($cstCB_flag) {
    push (@header, "CB-CB csts from ksync alignment");
} else {
    &abort ("must have CA-CA and/or CB/CB constraints");
}
push (@header, "measured in $pdb", $num_recs);

if ($cst_file) {
    open (OUTCST, '>'.$cst_file);
    select (OUTCST);
}
print join ("\n", @header, @constraints) . "\n";
if ($cst_file) {
    close (OUTCST);
    select (STDOUT);
}


# done
exit 0;

###############################################################################
# subs
###############################################################################

sub mapResCode {
    local $incode = shift;
    $incode = uc $incode;
    my $newcode = undef;

    my %one_to_three = ( 'G' => 'GLY',
                         'A' => 'ALA',
                         'V' => 'VAL',
                         'L' => 'LEU',
                         'I' => 'ILE',
                         'P' => 'PRO',
                         'C' => 'CYS',
                         'M' => 'MET',
                         'H' => 'HIS',
                         'F' => 'PHE',
                         'Y' => 'TYR',
                         'W' => 'TRP',
                         'N' => 'ASN',
                         'Q' => 'GLN',
                         'S' => 'SER',
                         'T' => 'THR',
                         'K' => 'LYS',
                         'R' => 'ARG',
                         'D' => 'ASP',
                         'E' => 'GLU',
                         'X' => 'XXX',
                         '0' => '  A',
                         '1' => '  C',
                         '2' => '  G',
                         '3' => '  T',
                         '4' => '  U'
                        );

    my %three_to_one = ( 'GLY' => 'G',
                         'ALA' => 'A',
                         'VAL' => 'V',
                         'LEU' => 'L',
                         'ILE' => 'I',
                         'PRO' => 'P',
                         'CYS' => 'C',
                         'MET' => 'M',
                         'HIS' => 'H',
                         'PHE' => 'F',
                         'TYR' => 'Y',
                         'TRP' => 'W',
                         'ASN' => 'N',
                         'GLN' => 'Q',
                         'SER' => 'S',
                         'THR' => 'T',
                         'LYS' => 'K',
                         'ARG' => 'R',
                         'ASP' => 'D',
                         'GLU' => 'E',

                         'UNK' => 'X',
                         '  X' => 'X',
                         '  A' => '0',
                         '  C' => '1',
                         '  G' => '2',
                         '  T' => '3',
                         '  U' => '4',
                         ' +A' => '0',
                         ' +C' => '1',
                         ' +G' => '2',
                         ' +T' => '3',
                         ' +U' => '4',
                        );


    my %fullname_to_one = ( 'GLYCINE'          => 'G',
                            'ALANINE'          => 'A',
                            'VALINE'           => 'V',
                            'LEUCINE'          => 'L',
                            'ISOLEUCINE'       => 'I',
                            'PROLINE'          => 'P',
                            'CYSTEINE'         => 'C',
                            'METHIONINE'       => 'M',
                            'HISTIDINE'        => 'H',
                            'PHENYLALANINE'    => 'F',
                            'TYROSINE'         => 'Y',
                            'TRYPTOPHAN'       => 'W',
                            'ASPARAGINE'       => 'N',
                            'GLUTAMINE'        => 'Q',
                            'SERINE'           => 'S',
                            'THREONINE'        => 'T',
                            'LYSINE'           => 'K',
                            'ARGININE'         => 'R',
                            'ASPARTATE'        => 'D',
                            'GLUTAMATE'        => 'E',
                            'ASPARTIC ACID'    => 'D',
                            'GLUTAMATIC ACID'  => 'E',
                            'ASPARTIC_ACID'    => 'D',
                            'GLUTAMATIC_ACID'  => 'E',
                            'SELENOMETHIONINE' => 'M',
                            'SELENOCYSTEINE'   => 'M',
                            'ADENINE'          => '0',
                            'CYTOSINE'         => '1',
                            'GUANINE'          => '2',
                            'THYMINE'          => '3',
                            'URACIL'           => '4'
                          );

    # map it
    #
    if (length $incode == 1) {
        $newcode = $one_to_three{$incode};
    }
    elsif (length $incode == 3) {
        $newcode = $three_to_one{$incode};
    }
    else {
        $newcode = $fullname_to_one{$incode};
    }


    # check for weirdness
    #
    if (! defined $newcode) {
#	&abort ("unknown residue '$incode'");
        print STDERR ("unknown residue '$incode' (mapping to 'Z')\n");
        $newcode = 'Z';
    }
    elsif ($newcode eq 'X') {
        print STDERR ("strange residue '$incode' (seen code, mapping to 'X')\n");
    }

    return $newcode;
}                    

# getCbCoords()
#
sub getCbCoords {
    my ($N_coords, $Ca_coords, $C_coords, $restype) = @_;
    my $Cb_coords = [];

    # formula (note: all vectors are normalized)
    # CaCb = bondlen * [cos5475*(-CaN -CaC) + sin5475*(CaN x CaC)]

    # config
    my $cos5475 = 0.577145190;                        # cos 54.75 = cos 109.5/2
    my $sin5475 = 0.816641555;                        # sin 54.75 = sin 109.5/2
    my $CC_bond = 1.536;                                          # from ethane
    my %CaCb_bond = ( 'A' => 1.524,
		      'C' => 1.531,
		      'D' => 1.532,
		      'E' => 1.530,
		      'F' => 1.533,
		      'G' => 1.532,
		      'H' => 1.533,
		      'I' => 1.547,
		      'K' => 1.530,
		      'L' => 1.532,
		      'M' => 1.530,
		      'N' => 1.532,
		      'P' => 1.528,
		      'Q' => 1.530,
		      'R' => 1.530,
		      'S' => 1.530,
		      'T' => 1.545,
		      'V' => 1.546,
		      'W' => 1.533,
		      'Y' => 1.533,
                      'ALA' => 1.524,
		      'CYS' => 1.531,
		      'ASP' => 1.532,
		      'GLU' => 1.530,
		      'PHE' => 1.533,
		      'GLY' => 1.532,
		      'HIS' => 1.533,
		      'ILE' => 1.547,
		      'LYS' => 1.530,
		      'LEU' => 1.532,
		      'MET' => 1.530,
		      'ASN' => 1.532,
		      'PRO' => 1.528,
		      'GLN' => 1.530,
		      'ARG' => 1.530,
		      'SER' => 1.530,
		      'THR' => 1.545,
		      'VAL' => 1.546,
		      'TRP' => 1.533,
		      'TYR' => 1.533,
		    );
    my $bondlen = (defined $restype) ? $CaCb_bond{$restype} : $CC_bond;

    # init vectors
    my $CaN  = +[];  my $CaN_mag  = 0.0;
    my $CaC  = +[];  my $CaC_mag  = 0.0;
    my $vert = +[];  my $vert_mag = 0.0;
    my $perp = +[];  my $perp_mag = 0.0;
    my $CaCb = +[];

    # CaN
    for ($i=0; $i<3; ++$i) {
	$CaN->[$i] = $N_coords->[$i] - $Ca_coords->[$i];
	$CaN_mag += $CaN->[$i] * $CaN->[$i];
    }
    $CaN_mag = sqrt ($CaN_mag);
    for ($i=0; $i<3; ++$i) {
	$CaN->[$i] /= $CaN_mag;
    }

    # CaC
    for ($i=0; $i<3; ++$i) {
	$CaC->[$i] = $C_coords->[$i] - $Ca_coords->[$i];
	$CaC_mag += $CaC->[$i] * $CaC->[$i];
    }
    $CaC_mag = sqrt ($CaC_mag);
    for ($i=0; $i<3; ++$i) {
	$CaC->[$i] /= $CaC_mag;
    }

    # vert = -CaN -CaC
    for ($i=0; $i<3; ++$i) {
	$vert->[$i] = - $CaN->[$i] - $CaC->[$i];
	$vert_mag += $vert->[$i] * $vert->[$i];
    }
    $vert_mag = sqrt ($vert_mag);
    for ($i=0; $i<3; ++$i) {
	$vert->[$i] /= $vert_mag;
    }

    # perp = CaN x CaC
    $perp->[0] = $CaN->[1] * $CaC->[2] - $CaN->[2] * $CaC->[1];
    $perp->[1] = $CaN->[2] * $CaC->[0] - $CaN->[0] * $CaC->[2];
    $perp->[2] = $CaN->[0] * $CaC->[1] - $CaN->[1] * $CaC->[0];
    # x product of two unit vectors is already unit, so no need to normalize

    # CaCb
    for ($i=0; $i<3; ++$i) {
	$CaCb->[$i] = $bondlen * ($cos5475 * $vert->[$i] +
				  $sin5475 * $perp->[$i]);
    }

    # Cb_coords
    for ($i=0; $i<3; ++$i) {
	$Cb_coords->[$i] = $Ca_coords->[$i] + $CaCb->[$i];
    }

    return $Cb_coords;
}

# typeI()
#
sub typeI {
    my $atomtype = shift;
    $atomtype =~ s/\s+//g;

    my %typenum = ( 'N'  => 0,
                    'CA' => 1,
                    'C'  => 2,
                    'O'  => 3,
                    'CB' => 4
                  );

    return $typenum{$atomtype};
}

# measureDistSq()
#
sub measureDistSq {
    my ($res_1, $res_2, $atomtype_1, $atomtype_2) = @_;
    my $atomtype_1_i = &typeI($atomtype_1);
    my $atomtype_2_i = &typeI($atomtype_2);
    
    &abort ("attempt to measure non-existent atoms")  if (! defined $res_1->[$atomtype_1_i] || ! defined $res_2->[$atomtype_2_i]);

    my $x_diff = $res_1->[$atomtype_1_i]->[0] - $res_2->[$atomtype_2_i]->[0];
    my $y_diff = $res_1->[$atomtype_1_i]->[1] - $res_2->[$atomtype_2_i]->[1];
    my $z_diff = $res_1->[$atomtype_1_i]->[2] - $res_2->[$atomtype_2_i]->[2];

    return ($x_diff*$x_diff + $y_diff*$y_diff + $z_diff*$z_diff);
}

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{usage: $0
\t -pdb           <pdb>
\t[-cstout        <cst_fileout>]     (def: STDOUT)
\t[-radius        <cst_radius>]      (def: 10.0 angstroms)
\t[-pad           <cst_pad>]         (def: 2.0 angstroms)
\t[-isep          <cst_isep>]        (def: 5 residues)
\t[-cstCA         <T/F>]             (def: F)
\t[-cstCB         <T/F>]             (def: T)
\t[-addlongrange  <T/F>]             (def: T)
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "pdb=s", 
		 "cstout=s", 
		 "radius=f", 
		 "pad=f", 
		 "isep=s", 
		 "cstCA=s", 
		 "cstCB=s", 
		 "addlongrange=s"
		);

    # Check for legal invocation
    if (! defined $opts{pdb}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{pdb});

    # Defaults
    $opts{radius}       = $def_cst_radius  if (! defined $opts{radius});
    $opts{pad}          = $def_cst_pad     if (! defined $opts{pad});
    $opts{isep}         = $def_cst_isep    if (! defined $opts{isep});
    #$opts{cstCA}        = 'TRUE'           if (! defined $opts{cstCA});
    $opts{cstCA}        = undef            if ($opts{cstCA} =~ /^F/i);
    $opts{cstCB}        = 'TRUE'           if (! defined $opts{cstCB});
    $opts{cstCB}        = undef            if ($opts{cstCB} =~ /^F/i);
    $opts{addlongrange} = 'TRUE'           if (! defined $opts{addlongrange});
    $opts{addlongrange} = undef            if ($opts{addlongrange} =~ /^F/i);

    return %opts;
}

###############################################################################
# util
###############################################################################

# insertSortIndexList ()
#
sub insertSortIndexList {
    my ($val_list, $direction) = @_;
    my $index_list = +[];
    my ($index, $val, $i, $i2, $assigned);

    $index_list->[0] = 0;
    for ($index=1; $index <= $#{$val_list}; ++$index) {
        $assigned = undef;
        $val = $val_list->[$index];
        for ($i=0; $i <= $#{$index_list}; ++$i) {
            if ($direction eq 'decreasing') {
                if ($val > $val_list->[$index_list->[$i]]) {
                    for ($i2=$#{$index_list}; $i2 >= $i; --$i2) {
                        $index_list->[$i2+1] = $index_list->[$i2];
                    }
                    $index_list->[$i] = $index;
                    $assigned = 'true';
                    last;
                }
            }
            else {
                if ($val < $val_list->[$index_list->[$i]]) {
                    for ($i2=$#{$index_list}; $i2 >= $i; --$i2) {
                        $index_list->[$i2+1] = $index_list->[$i2];
                    }
                    $index_list->[$i] = $index;
                    $assigned = 'true';
                    last;
                }
            }
        }
        $index_list->[$#{$index_list}+1] = $index  if (! $assigned);
    }
    return $index_list;
}
     
# runCmd
#
sub runCmd {
    my $cmd = shift;
    #print $cmd."\n" if ($debug);
    &abort ("failure running cmd: $cmd\n") if (system ($cmd) != 0);
    #if (system ($cmd) != 0) {
    #    print STDERR ("failure running cmd: $cmd\n");
    #    return -2;
    #}
    return 0;
}

# logMsg()
#
sub logMsg {
    my ($msg, $logfile) = @_;

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    else {
	select (STDERR);
    }
    print $msg, "\n";
    if ($logfile) {
        close (LOGFILE);
    }
    select (STDOUT);

    return 'true';
}

# checkExist()
#
sub checkExist {
    my ($type, $path) = @_;
    if ($type eq 'd') {
	if (! -d $path) { 
            print STDERR "$0: dirnotfound: $path\n";
            exit -3;
	}
    }
    elsif ($type eq 'f') {
	if (! -f $path) {
            print STDERR "$0: filenotfound: $path\n";
            exit -3;
	}
	elsif (! -s $path) {
            print STDERR "$0: emptyfile: $path\n";
            exit -3;
	}
    }
}

# abort()
#
sub abort {
    my $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}

# writeBufToFile()
#
sub writeBufToFile {
    my ($file, $bufptr) = @_;
    if (! open (FILE, '>'.$file)) {
	&abort ("$0: unable to open file $file for writing");
    }
    print FILE join ("\n", @{$bufptr}), "\n";
    close (FILE);
    return;
}

# fileBufString()
#
sub fileBufString {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

# fileBufArray()
#
sub fileBufArray {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

# bigFileBufArray()
#
sub bigFileBufArray {
    my $file = shift;
    my $buf = +[];
    if ($file =~ /\.gz|\.Z/) {
        if (! open (FILE, "gzip -dc $file |")) {
            &abort ("$0: unable to open file $file for gzip -dc");
        }
    }
    elsif (! open (FILE, $file)) {
        &abort ("$0: unable to open file $file for reading");
    }
    while (<FILE>) {
        chomp;
        push (@$buf, $_);
    }
    close (FILE);
    return $buf;
}     

###############################################################################
# end
1;                                                     # in case it's a package
###############################################################################
