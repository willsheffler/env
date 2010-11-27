#!/usr/bin/perl
##
## Copyright 2003, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.4 $
##  $Date: 2004/07/15 23:04:57 $
##  $Author: dylan $
##
###############################################################################

# timing trick
#
#$timing          = 0.0;
#$timing_interval = 0.100;  # 100 millisecs
#$max_timing      = 3.0;
#if (fork() == 0) {
#    # then we're child
#    while ($timing < $max_timing) {
#	print "TIME: $timing\n";
#	select(undef, undef, undef, $timing_interval); 
#	$timing += $timing_interval;
#    }
#    exit 0;
#}
			       
###############################################################################
# conf
###############################################################################

# general
#
$| = 1;                                              # disable stdout buffering
$debug = 0;                                          # 1: chatter while running
$use_filesys = undef;                     # don't output mammoth to file system
$use_filesys = 'true';  # DEBUG                # output mammoth to file system
#$use_tmp_pdbs = undef;

# paths
if (! defined $ENV{'BAKER_HOME'}) {
    &abort ("must define \$BAKER_HOME in environment");
}
if (! defined $ENV{'SCRATCH_DIR'}) {
    &abort ("must define \$SCRATCH_DIR in environment");
}
$src_dir = $ENV{'BAKER_HOME'}."/src";
$dat_dir = $ENV{'BAKER_HOME'}."/dat";
$tmp_dir = (-d $ENV{'SCRATCH_DIR'}."/tmp") ? $ENV{'SCRATCH_DIR'}."/tmp" : "/tmp";   

$mammoth = "$src_dir/shareware/mastodon/mammoth";
&abort ("can't find program $mammoth")  if (! -f $mammoth);


# vars
#
@gdt_devs  = (1.00, 2.0, 4.0,  8.0);
@gdt_rmsds = (0.75, 1.5, 3.0,  6.0);         # big enough not to bound (r=3d/4)

@grt_rmsds = (1.00, 2.0, 4.0,  6.0);
@grt_devs  = (1.75, 3.5, 7.0, 10.5);    # not too big: kills maxsub (d=2r-2r/8)


###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
my $native_pdb_file = $opts{nativepdbfile};
my $model_pdb_file  = $opts{modelpdbfile};
my $out_file        = $opts{outfile};
my $dev_thresh      = $opts{devthresh};
my $rmsd_thresh     = $opts{rmsdthresh};
my $verbose         = $opts{verbose};
my $wrt_common      = $opts{common};
my $wrt_both        = $opts{both};
my $perc_flag       = $opts{perc};
my $gdt_mm          = $opts{gdtmm};
my $grt_mm          = $opts{grtmm};
my $seq_independent = $opts{seqindependent};

if ($native_pdb_file =~ /\.gz$|\.Z$/ || $model_pdb_file =~ /\.gz$|\.Z$/) {
    $use_tmp_pdbs = 'true';
}

my $basename_native = $native_pdb_file;
$basename_native    =~ s!^.*/!!;
$basename_native    =~ s!\..*$!!;
my $basename_model  = $model_pdb_file;
$basename_model     =~ s!^.*/!!;
$basename_model     =~ s!\..*$!!;

###############################################################################
# main
###############################################################################

# clean native pdb
#
@native_res             = ();
@native_res_n           = ();
@native_denseCA_pdb_buf = ();
@native_pdb_buf = &fileBufArray ($native_pdb_file);
for ($i=0; $i <= $#native_pdb_buf; ++$i) {
    $line = $native_pdb_buf[$i];
    if ($line =~ /^ATOM/) {
        next if (substr ($line, 12, 4) ne " CA ");
        $occ = substr ($line, 54, 6);
        next if ($occ < 0.0);
	$res_n = substr ($line, 22, 4);
	$res_n =~ s/\s+//g;

        substr ($line, 21, 1) = ' ';                            # remove chain
        $native_denseCA_pdb_buf[$res_n] = $line;

	$native_res[$res_n] = substr ($line, 17, 3);
	push (@native_res_n, $res_n);
	++$native_len;
    }
}


# clean model pdb
#
@model_res             = ();
@model_res_n           = ();
@model_denseCA_pdb_buf = ();
@model_pdb_buf = &fileBufArray ($model_pdb_file);
for ($i=0; $i <= $#model_pdb_buf; ++$i) {
    $line = $model_pdb_buf[$i];
    if ($line =~ /^ATOM/) {
        next if (substr ($line, 12, 4) ne " CA ");
        $occ = substr ($line, 54, 6);
        next if ($occ < 0.0);
	$res_n = substr ($line, 22, 4);
	$res_n =~ s/\s+//g;

        substr ($line, 21, 1) = ' ';                            # remove chain
        $model_denseCA_pdb_buf[$res_n] = $line;

	$model_res[$res_n] = substr ($line, 17, 3);
	push (@model_res_n, $res_n);
    }
}


# determine native length and common length
#
$native_len = 0;
$common_len = 0;
for ($res_n=1; $res_n <= $#native_res; ++$res_n) {
    ++$native_len  if (defined $native_res[$res_n]);
    if (defined $native_res[$res_n] && defined $model_res[$res_n]) {
	&abort ("non-equivalent residues at residue number $res_n: NATIVE: $native_res[$res_n] MODEL: $model_res[$res_n]")  if ($native_res[$res_n] ne $model_res[$res_n]);

	++$common_len;
    }
}
&abort ("no positions in common")  if ($common_len == 0);


# get multiple scores if GDT_MM or GRT_MM
#
@rmsds = ();
@devs  = ();
if ($gdt_mm) {
    @rmsds = (@gdt_rmsds, $rmsd_thresh);
    @devs  = (@gdt_devs,  $dev_thresh);
}
elsif ($grt_mm) {
    @rmsds = (@grt_rmsds, $rmsd_thresh);
    @devs  = (@grt_devs,  $dev_thresh);
}
else {
    @rmsds = ($rmsd_thresh);
    @devs  = ($dev_thresh);
}


# if our pdbs are zipped, then we'll operate on local copies
#
if ($use_tmp_pdbs) {
    $native_tmp_pdb_file = "$tmp_dir/$basename_native";
    open (TMP_PDB, '>'.$native_tmp_pdb_file);
    print TMP_PDB join ("\n", @native_denseCA_pdb_buf)."\n";
    close (TMP_PDB);
    $model_tmp_pdb_file = "$tmp_dir/$basename_model";
    open (TMP_PDB, '>'.$model_tmp_pdb_file);
    print TMP_PDB join ("\n", @model_denseCA_pdb_buf)."\n";
    close (TMP_PDB);
    $use_native_pdb_file = $native_tmp_pdb_file;
    $use_model_pdb_file  = $model_tmp_pdb_file;
} else {
    $use_native_pdb_file = $native_pdb_file;
    $use_model_pdb_file  = $model_pdb_file;
}


# now get Mammoth MaxsubRMSDs at each threshold (last one is requested)
#
@range_scores = ();
for ($thresh_i=0; $thresh_i <= $#rmsds; ++$thresh_i) {

    # mammoth align native and model
    #
    @mammoth_out_buf = ();
    $just_maxsub_opt = ($seq_independent) ? "" : "-M 1";
    if ($use_filesys) {
	$mammoth_out_file = "$tmp_dir/$basename_native-$basename_model.mammoth_maxsubRMSD-$rmsds[$thresh_i]-$devs[$thresh_i]";
	&runCmd (qq{$mammoth -p $use_model_pdb_file -e $use_native_pdb_file -o $mammoth_out_file -R $rmsds[$thresh_i] -D $devs[$thresh_i] $just_maxsub_opt -r 0 2> /dev/null});
	@mammoth_out_buf = &fileBufArray ($mammoth_out_file);
    } else {
	@mammoth_out_buf = split (/\n/, `$mammoth -p $use_model_pdb_file -e $use_native_pdb_file -R $rmsds[$thresh_i] -D $devs[$thresh_i] $just_maxsub_opt -r 0 2> /dev/null`);
    } 
    unlink "rasmol.tcl"       if (-f "rasmol.tcl");       # shouldn't need: -r0
    unlink "maxsub_sup.pdb"   if (-f "maxsub_sup.pdb");   # shouldn't need: -r0
    unlink "maxsub_sup2.pdb"  if (-f "maxsub_sup2.pdb");  # shouldn't need: -r0


    # extract maxsub scores  
    #
    $maxsubRMSD_wrt_native = undef;
    $maxsubRMSD_wrt_common = undef;
    $maxsubRMSD_wrt_both   = undef;
    $nali                  = undef;
    $rmsd                  = undef;
    foreach $line (@mammoth_out_buf) {
	if ($line =~ /^\s*PSI\(end\)\s*=\s*([\d\.]+)\s*NALI\s*=\s*(\d+)\s*NORM\s*=\s*(\d+)\s*RMS\s*=\s*([\d\.]+)/) {
	    $psi_score = $1;
	    $nali      = $2;
	    $norm      = $3;
	    $rmsd      = $4;
	    
	    # don't trust nali because mammoth often/always? drops last residue)
	    #$maxsubRMSD_wrt_native = $nali / $native_len;
	    #$maxsubRMSD_wrt_common = $nali / $common_len;
	    #$maxsubRMSD_wrt_both   = ($maxsubRMSD_wrt_native+$maxsubRMSD_wrt_common) / 2.0;
	    
	    last;
	}
    }


    # read alignment, and calculate nali ourselves
    #
    $model_res_i  = -1;
    $native_res_i = -1;
    $last_aligned_model_res_i  = -1;
    $last_aligned_native_res_i = -1;
    $nali                =  0;
    @align_map           = ();
    for ($line_i=0; $line_i <= $#mammoth_out_buf; ++$line_i) {
	if ($mammoth_out_buf[$line_i] =~ /^Prediction /) {
	    $aligned_str   = substr ($mammoth_out_buf[$line_i-1], 11), 
	    $model_aln_str = substr ($mammoth_out_buf[$line_i],   11);
	    ++$line_i;
	}
	elsif ($mammoth_out_buf[$line_i] =~ /^Experiment /) {
	    $native_aln_str = substr ($mammoth_out_buf[$line_i+1], 11);
	    ++$line_i;
	    
	    @aligned    = split (//, $aligned_str);
	    @model_aln  = split (//, $model_aln_str);
	    @native_aln = split (//, $native_aln_str);
	    
	    for ($pos_i=0; $pos_i <= $#model_aln; ++$pos_i) {
		if ($model_aln[$pos_i] ne ' ') {
		    ++$model_res_i   if ($model_aln[$pos_i]  ne '.');
		    ++$native_res_i  if ($native_aln[$pos_i] ne '.');
		    if ($aligned[$pos_i] eq '*') {
			if ($seq_independent ||
			    $native_res_n[$native_res_i] == $model_res_n[$model_res_i]) {
			    $align_map[$model_res_i] = $native_res_i;
			    ++$nali;
			    
			    $last_aligned_model_res_i  = $model_res_i;
			    $last_aligned_native_res_i = $native_res_i;
			}
		    }
		}
	    }
	}
    }
    # if next to last residue is last reported by mammoth, and is aligned,
    #   then assume non-reported but real last residue is also aligned
    if (($last_aligned_model_res_i  == $#model_res_n-1 && $last_aligned_model_res_i == $model_res_i) ||
	($last_aligned_native_res_i == $#native_res_n-1 && $last_aligned_native_res_i == $native_res_i)) {
	
	if ($seq_independent ||
	    $native_res_n[$last_aligned_native_res_i] == $model_res_n[$last_aligned_model_res_i]) {
	    
	    ++$nali;
	    ++$last_aligned_model_res_i;
	    ++$last_aligned_native_res_i;
	    $align_map[$last_aligned_model_res_i] = $last_aligned_native_res_i;
	}
    }
    
    
    # determine zones
    #
    if ($thresh_i == $#rmsds) {
	@zones = ();
	for ($model_res_i=0, $start_model_res_i=undef, $last_native_res_i=-1000;
	     $model_res_i <= $#align_map+1;
	     ++$model_res_i) {
	    
	    $native_res_i = $align_map[$model_res_i];
	    if (! defined $native_res_i || 
		(defined $native_res_i && $native_res_i != $last_native_res_i+1)
		) {
		if (defined $start_model_res_i) {
		    $start_model_res_n  = $model_res_n[$start_model_res_i];
		    $last_model_res_n   = $model_res_n[$last_model_res_i];
		    $start_native_res_n = $native_res_n[$start_native_res_i]; 
		    $last_native_res_n  = $native_res_n[$last_native_res_i];
		    
		    push (@zones, sprintf (qq{zone %4d-%-4d:%4d-%-4d}, 
					   $start_model_res_n, 
					   $last_model_res_n, 
					   $start_native_res_n, 
					   $last_native_res_n));
		}
		if (defined $native_res_i) {
		    $start_model_res_i  = $model_res_i;
		    $start_native_res_i = $native_res_i;
		} else {
		    $start_model_res_i  = undef;
		    $start_native_res_i = undef;
		}
	    }
	    $last_model_res_i = $model_res_i;
	    $last_native_res_i = $native_res_i  if (defined $native_res_i);
	}                
    }              
    
    
    # now calculate maxsubRMSD
    #
    $maxsubRMSD_wrt_native = $nali / $native_len;
    $maxsubRMSD_wrt_common = $nali / $common_len;
    $maxsubRMSD_wrt_both   = ($maxsubRMSD_wrt_native+$maxsubRMSD_wrt_common) / 2.0;

    if (($gdt_mm || $grt_mm) && $thresh_i != $#rmsds) {
	if ($wrt_common) {
	    push (@range_scores, $maxsubRMSD_wrt_common);
	} 
	elsif ($wrt_both) {
	    push (@range_scores, $maxsubRMSD_wrt_both);
	}
	else {
	    push (@range_scores, $maxsubRMSD_wrt_native);
	}
    }


    # clean up
    #
    unlink $mammoth_out  if (-f $mammoth_out);
}


# clean up
#
if ($use_tmp_pdbs) {
    unlink $native_tmp_pdb_file;
    unlink $model_tmp_pdb_file;
}


# get GDT_MM or GRT_MM scores
#
if ($gdt_mm || $grt_mm) {
    $range_score = 0.0;
    foreach $score (@range_scores) {
	$range_score += $score;
    }
    $range_score /= 4.0;

    if ($gdt_mm) {
	$gdt_mm_score = sprintf ("%-6.4f", $range_score);
	$grt_mm_score = 'NA';
    }
    else {
	$gdt_mm_score = 'NA';
	$grt_mm_score = sprintf ("%-6.4f", $range_score);
    }
} else {
    $gdt_mm_score = 'NA';
    $grt_mm_score = 'NA';
}


# report scores
#
&abort ("couldn't find score line")  if (! defined $rmsd);
@out_buf = ();
$native_pdb_name = $native_pdb_file;
$native_pdb_name =~ s!^.*\/!!;
$model_pdb_name  = $model_pdb_file;
$model_pdb_name  =~ s!^.*\/!!;

if ($verbose) {
    push (@out_buf, sprintf ("native:                 %-s", $native_pdb_name));
    push (@out_buf, sprintf ("model:                  %-s", $model_pdb_name));
    push (@out_buf, sprintf ("dev_thresh:             %-6.2f", $dev_thresh));
    push (@out_buf, sprintf ("rmsd_thresh:            %-6.2f", $rmsd_thresh));
    push (@out_buf, sprintf ("native_len:             %-4d",   $native_len));
    push (@out_buf, sprintf ("common_len:             %-4d",   $common_len));
    push (@out_buf, sprintf ("fit_len:                %-4d",   $nali));
    push (@out_buf, sprintf ("fit_rmsd:               %-6.2f", $rmsd));
    push (@out_buf, sprintf ("maxsubRMSD_wrt_native:  %-6.4f", $maxsubRMSD_wrt_native));
    push (@out_buf, sprintf ("maxsubRMSD_wrt_common:  %-6.4f", $maxsubRMSD_wrt_common));
    push (@out_buf, sprintf ("maxsubRMSD_wrt_both:    %-6.4f", $maxsubRMSD_wrt_both));
    push (@out_buf, sprintf ("GDT_MM:                 %-6s", $gdt_mm_score));
    push (@out_buf, sprintf ("GRT_MM:                 %-6s", $grt_mm_score));
    push (@out_buf, sprintf ("REFERENCE:              %-s", (($wrt_common) ? "COMMON" : (($wrt_both) ? "BOTH" : "NATIVE"))));
    push (@out_buf, sprintf ("SEQ-(IN)DEPENDENT:      %-s", (($seq_independent) ? 'SEQ-INDEP' : 'SEQ-DEP')));
    push (@out_buf, "ZONES:");
    push (@out_buf, @zones);
} else {
    push (@out_buf, join ("\t", 
			  $native_pdb_name,
			  $model_pdb_name,
			  sprintf ("%-4d", $nali),
			  sprintf ("%-6.2f", $rmsd),
			  sprintf ("%-6.4f", $maxsubRMSD_wrt_native),
			  sprintf ("%-6.4f", $maxsubRMSD_wrt_common),
			  sprintf ("%-6.4f", $maxsubRMSD_wrt_both),
			  sprintf ("%-6s", $gdt_mm_score),
			  sprintf ("%-6s", $grt_mm_score),
			  sprintf ("%-4d", $native_len),
			  sprintf ("%-4d", $common_len),
			  sprintf ("%-6.2f", $dev_thresh),
			  sprintf ("%-6.2f", $rmsd_thresh),
			  (($wrt_common) ? "COMMON" : (($wrt_both) ? "BOTH" : "NATIVE")),
			  (($seq_independent) ? 'SEQ-INDEP' : 'SEQ-DEP')
			 )
	  );
}
if ($out_file) {
    open (OUT, '>'.$out_file);
    select (OUT);
}
print join ("\n", @out_buf)."\n";
if ($out_file) {
    close (OUT);
    select (STDOUT);
}


# done
exit 0;

###############################################################################
# subs
###############################################################################

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{
usage: $0
\t -nativepdbfile  <native_pdb_file>
\t -modelpdbfile   <model_pdb_file>
\t[-outfile        <out_file>]       (def: STDOUT)
\t[-devthresh      <dev_thresh>]     (def: 7.0 A)
\t[-rmsdthresh     <rmsd_thresh>]    (def: 4.0 A)
\t[-verbose]                         (def: undef, i.e. brief)
\t[-common or -both]                 (def: undef, i.e. wrt native)
\t[-gdtmm or -grtmm]                 (def: undef, i.e. don\'t get GDT_MM or GRT_MM)
\t[-seqindependent]                  (def: undef, i.e. sequence dependent)

NON-VERBOSE OUTPUT LINE FORMAT (tab separated):
\tnative_pdb_file     (native)
\tmodel_pdb_file      (model)
\tn_ali               (raw count of aligned residues)
\trmsd_ali            (rmsd of aligned residues)
\tmaxsubRMSD_native   (percentage fit wrt native)
\tmaxsubRMSD_common   (percentage fit wrt model\'s common residues with native)
\tmaxsubRMSD_both     (average of native and common)
\tGDT_MM              (global distance test - mammoth maxsubRMSD: 1,2,4,8 A)
\tGRT_MM              (global rmsd test     - mammoth maxsubRMSD: 1,2,4,6 A)
\tnative_len          (number of occupied CA atoms in native)
\tcommon_len          (number of occupied CA atoms in model and native)
\tdev_thresh          (deviation threshold for maxsubRMSD)
\trmsd_thresh         (rmsd threshold for maxsubRMSD)
\tG[DR]T_MM_reference (GDT_MM or GRT_MM with respect to NATIVE/COMMON/BOTH)
\tseq_(in)dep         (sequence-dependent or sequence-independent alignment)
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "nativepdbfile=s", 
		 "modelpdbfile=s",
		 "outfile=s",
		 "devthresh=f",
		 "rmsdthresh=f",
		 "verbose",
		 "common",
		 "both",
                 "gdtmm",
                 "grtmm",
		 "seqindependent");

    # Check for legal invocation
    if (! defined $opts{nativepdbfile} ||
	! defined $opts{modelpdbfile}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{nativepdbfile});
    &checkExist ('f', $opts{modelpdbfile});

    # Complain about weird command lines
    &abort ("can't report wrt common and both") if (defined $opts{common} && defined $opts{both});
    &abort ("can't report both GDT_MM and GRT_MM") if (defined $opts{gdtmm} && defined $opts{grtmm});

    # Defaults
    $opts{devthresh}  = 7.0  if (! defined $opts{devthresh});
    $opts{rmsdthresh} = 4.0  if (! defined $opts{rmsdthresh});

    return %opts;
}

###############################################################################
# util
###############################################################################

# chompEnds ()
#
sub chompEnds {
    my $str = shift;
    $str =~ s/^\s+|\s+$//g;
    return $str;
}


# cleanStr ()
#
sub cleanStr {
    my $str = shift;
    $str =~ s/[\x00-\x08\x0B-\x1F\x80-\xFF]//g;
    return $str;
}


# iterElimSortIndexList ()
#
sub iterElimSortIndexList {
    my ($val1_list, $val2_list, $fraction, $direction) = @_;
    my $index_list        = +[];
    my $local_index_list  = +[];
    my $local_val_list    = +[];
    my $local_sorted_list = +[];
    my ($index, $i, $j);

    my $sorted_val1_list = &insertSortIndexList ($val1_list, $direction);
    for ($i=0; $i <= $#{$sorted_val1_list}; ++$i) {
	$index_list->[$i] = $sorted_val1_list->[$i];
    }

    my $done = undef;
    my $toggle = 2;
    $cut = int ($#{$index_list} * $fraction);
    $last_cut = $#{$index_list};
    while ($cut > 0) {
	# sort the right half ("discards")
	$local_index_list = +[];
	$local_val_list   = +[];
	for ($j=0; $cut+$j+1 <= $last_cut; ++$j) {
	    $index                  = $index_list->[$cut+$j+1];
	    $local_index_list->[$j] = $index;
	    $local_val_list->[$j]   = ($toggle == 1) ? $val1_list->[$index]
		                                     : $val2_list->[$index];
	}
	$local_sorted_index_list = &insertSortIndexList ($local_val_list, $direction);
	for ($j=0; $cut+$j+1 <= $last_cut; ++$j) {
	    $local_index = $local_sorted_index_list->[$j];
	    $index_list->[$cut+$j+1] = $local_index_list->[$local_index];
	}

	# sort the left half ("keeps")
	$local_index_list = +[];
	$local_val_list   = +[];
	for ($j=0; $j <= $cut; ++$j) {
	    $index                  = $index_list->[$j];
	    $local_index_list->[$j] = $index;
	    $local_val_list->[$j]   = ($toggle == 1) ? $val1_list->[$index]
		                                     : $val2_list->[$index];
	}
	$local_sorted_index_list = &insertSortIndexList ($local_val_list, $direction);
	for ($j=0; $j <= $cut; ++$j) {
	    $local_index = $local_sorted_index_list->[$j];
	    $index_list->[$j] = $local_index_list->[$local_index];
	}
	
	# update cut and toggle
	$toggle = ($toggle == 1) ? 2 : 1;
	$last_cut = $cut;
	$cut = int ($last_cut * $fraction);
    }

    return $index_list;
}

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

# readFiles
#
sub readFiles {
    my ($dir, $fullpath_flag) = @_;
    my $inode;
    my @inodes = ();
    my @files = ();
    
    opendir (DIR, $dir);
    @inodes = sort readdir (DIR);
    closedir (DIR);
    foreach $inode (@inodes) {
	next if (! -f "$dir/$inode");
	next if ($inode =~ /^\./);
	push (@files, ($fullpath_flag) ? "$dir/$inode" : "$inode");
    }
    return @files;
}

# createDir
#
sub createDir {
    my $dir = shift;
    if (! -d $dir && (system (qq{mkdir -p $dir}) != 0)) {
	print STDERR "$0: unable to mkdir -p $dir\n";
	exit -2;
    }
    return $dir;
}

# copyFile
#
sub copyFile {
    my ($src, $dst) = @_;
    if (-f $src) {
	if (system (qq{cp $src $dst}) != 0) {
	    print STDERR "$0: unable to cp $src $dst\n";
	    exit -2;
	}
    } else {
	print STDERR "$0: file not found: '$src'\n";
    }
    return $dst;
}

# zip
#
sub zip {
    my $file = shift;
    if ($file =~ /^\.Z/ || $file =~ /\.gz/) {
	&abort ("already a zipped file $file");
    }
    if (-s $file) {
	if (system (qq{gzip -9f $file}) != 0) {
	    &abort ("unable to gzip -9f $file");
	}
    } elsif (-f $file) {
	&abort ("file empty: '$file'");
    } else {
	&abort ("file not found: '$file'");
    }
    $file .= ".gz";
    return $file;
}

# unzip
#
sub unzip {
    my $file = shift;
    if ($file !~ /^\.Z/ && $file !~ /\.gz/) {
	&abort ("not a zipped file $file");
    }
    if (-f $file) {
	if (system (qq{gzip -d $file}) != 0) {
	    &abort ("unable to gzip -d $file");
	}
    } else {
	&abort ("file not found: '$file'");
    }
    $file =~ s/\.Z$|\.gz$//;
    if (! -s $file) {
	&abort ("file empty: '$file'");
    }
    return $file;
}

# remove
#
sub remove {
    my $inode = shift;
    if (-e $inode) {
	if (system (qq{rm -rf $inode}) != 0) {
	    print STDERR "$0: unable to rm -rf $inode\n";
	    exit -2;
	}
    } else {
	print STDERR "$0: inode not found: '$inode'\n";
    }
    return $inode;
}
     
# runCmd
#
sub runCmd {
    my ($cmd, $nodie, $silent) = @_;
    my $ret;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print "[$date]:$0:RUNNING: $cmd\n" if ($debug && ! $silent);
    $ret = system ($cmd);
    #$ret = ($?>>8)-256;
    if ($ret != 0) {
	$date = `date +'%Y-%m-%d_%T'`;  chomp $date;
	print STDERR ("[$date]:$0: FAILURE (exit: $ret): $cmd\n");
	if ($nodie) {
	    return $ret;
	} else {
	    exit $ret;
	}
    }
    return 0;
}

# logMsg()
#
sub logMsg {
    my ($msg, $logfile) = @_;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    else {
	select (STDERR);
    }
    print "[$date]:$0: $msg\n";
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

# alert()
#
sub alert {
    my $msg = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print STDERR "[$date][ALERT][$0] $msg\n";
    return;
}

# abort()
#
sub abort {
    my $msg = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print STDERR "[$date][ABORT][$0] $msg\n";
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
    if ($file =~ /\.gz$|\.Z$/) {
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
    if ($file =~ /\.gz$|\.Z$/) {
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
    if ($file =~ /\.gz$|\.Z$/) {
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
