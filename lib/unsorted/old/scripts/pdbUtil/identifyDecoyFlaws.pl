#!/usr/bin/perl
##
## Copyright 2003, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.11 $
##  $Date: 2004/07/15 23:04:57 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

# general
$| = 1;                                              # disable stdout buffering
$debug = 1;                                             # chatter while running

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

# vars
$CA_max_dist_per_res = 4.0;  # in angstroms (extended chain is more like 3.7)
$CA_CA_min_dist      = 4.0;  # in angstroms (extended chain is more like 3.7)

@target_fields = (qw(Eres Erep Eaa));
#@target_fields = (qw(Eaa));

%thresholds = ( 'Eres' => 5.0,
		'Erep' => 4.0,
		'Eaa'  => 1.5,
	      );

%flaw_codes = ( 'break' => 1,
		'crush' => 2,
		'Eres'  => 3,
		'Erep'  => 4,
		'Eaa'   => 5,
	      );

###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
my $zones_file     = $opts{zonesfile};
my $ssa_file       = $opts{ssafile};
my $complete_file  = $opts{completefile};
my $scored_file    = $opts{scoredfile};
my $out_file       = $opts{outfile};

#if ($#ARGV < 0) {
#    print STDERR "usage: $0 <file>\n";
#    exit -1;
#}
#$file = shift @ARGV;

$zones_buf    = &bigFileBufArray ($zones_file);
$complete_buf = &bigFileBufArray ($complete_file);
$scored_buf   = &bigFileBufArray ($scored_file);

###############################################################################
# main
###############################################################################


# read complete model ssa
#
@ssa = ();
foreach $line (&fileBufArray ($ssa_file)) {
    next if ($line =~ /^\s*>/);
    $line =~ s/\s+//g;
    push (@ssa, split (//, $line));
}


# read zones
#
@template_mask = ();
@flawed_gap    = ();
$last_q_stop   = undef;
$last_p_stop   = undef;
foreach $line (@{$zones_buf}) {
    next if ($line !~ /^\s*zone\s*(\d+)\s*-\s*(\d+)\s*:\s*(\d+)\s*-\s*(\d+)/);
    $q_start = $1;
    $q_stop  = $2;
    $p_start = $3;
    $p_stop  = $4;

    # simple check for flawed gap
    #
    if (defined $last_q_stop && $q_start == $last_q_stop+1 &&
	defined $last_p_stop && $p_start != $last_p_stop+1
	) {
	$flawed_gap[$q_start-1] = $flaw_codes{'break'};
    }
    if (defined $last_p_stop && $p_start == $last_p_stop+1 &&
	defined $last_q_stop && $q_start != $last_q_stop+1
	) {
	$flawed_gap[$q_start-1] = $flaw_codes{'crush'};
    }
    $last_q_stop = $q_stop;
    $last_p_stop = $p_stop;


    # assign template mask
    #
    for ($res_j=$q_start-1; $res_j <= $q_stop-1; ++$res_j) {
	$template_mask[$res_j] = 1;
    }
}


# check for reasonableness of physical breaks longer than zero residues, 
#   and store density_map
#
$started = undef;
$res_i = -1;
$res_j = -1;
@occ = ();
@density_map = ();
$last_occ_res_j = undef;
$coords = +[];
@query_fasta = ();
foreach $line (@$complete_buf) {
    if ($line =~ /^ATOM/) {
	$started = 'true';
    } elsif ($started) {
	last;
    }

    if (substr ($line, 12, 4) eq ' CA ') {
	++$res_j;
	$occ[$res_j] = substr ($line, 54, 6);
	$query_fasta[$res_j] = &mapResCode (substr ($line, 17, 3));
	if ($occ[$res_j] >= 0) {
	    $density_map[++$res_i] = $res_j;
	    $coords->[$res_j]->[0] = substr ($line, 30, 8);
	    $coords->[$res_j]->[1] = substr ($line, 38, 8);
	    $coords->[$res_j]->[2] = substr ($line, 46, 8);

	    if (defined $last_occ_res_j) {
		$q_stem_start = $last_occ_res_j;
		$q_stem_stop  = $res_j;
		$loop_len = $q_stem_stop - $q_stem_start;
		$max_dist = ($loop_len+1) * $CA_max_dist_per_res;
		if (&coord_dist_sq ($coords->[$res_j], $coords->[$last_occ_res_j]) > $max_dist*$max_dist) {
		    $flawed_gap[$res_j] = $flaw_codes{'break'};
		}
	    }
	    $last_occ_res_j = $res_j;
	}
    }
}


# read per residue energies
#
$started = undef;
$started_avg = undef;
%field_col = ();
@energies = ();
$res_i = -1;
@flawed_res = ();
foreach $line (@$scored_buf) {
    if ($avg && $line =~ /^\s*energies-average/) {
	$started_avg = 'true';
	next;
    }
    next if ($avg && ! $started_avg);

    if ($line =~ /^\s*res/) {
	$line =~ s/^\s+|\s+$//g;
	@fields = split (/\s+/, $line);
	for ($header_i=0; $header_i <= $#fields; ++$header_i) {
	    $field_col{$fields[$header_i]} = $header_i;
	}
	$started = 'true';
	next;
    }

    next if (! $started);
    last if ($line !~ /^\s*\d+\s*/);

    ++$res_i;
    $res_j = $density_map[$res_i];

    $line =~ s/^\s+|\s+$//g;
    @scores = split (/\s+/, $line);
    for ($field_i=0; $field_i <= $#target_fields; ++$field_i) {
	$field = $target_fields[$field_i];
	$energy = $scores[$field_col{$field}];
	$energies[$res_i]->{$field} = $energy;
	$flawed_res[$res_j] = $flaw_codes{$field}  if ($energy >= $thresholds{$field});
# DEBUG
#	print "$res_i $field $energy\n";
#	if ($energy >= $thresholds{$field}) {
#	    print "\tSTORING\n";
#	}
    }
}


# build out_buf
#
@out_buf = ();
for ($res_j=0; $res_j <= $#query_fasta; ++$res_j) {
    $out_buf[0] .= $query_fasta[$res_j];
    $out_buf[2] .= ($template_mask[$res_j]) ? 'X' : '.';
    if ($occ[$res_j] >= 0) {
	$out_buf[1] .= $ssa[$res_j];
	$flaw = (defined $flawed_gap[$res_j]) ? $flawed_gap[$res_j] : $flawed_res[$res_j];	
	$out_buf[3] .= (defined $flaw) ? $flaw : '-';
#	print STDERR "flawed res $query_fasta[$res_j]".($res_j+1). " flaw: $flaw\n"  if (defined $flaw);

    } else {
	$out_buf[1] .= '.';
	$out_buf[3] .= '.';
    }
}


# output
#
if ($out_file) {
    open (OUT, '>'.$out_file);
    select (OUT);
}
print join ("\n", @out_buf) ."\n"; 
if ($out_file) {
    close (OUT);
    select (STDOUT);
}


# done
exit 0;

###############################################################################
# subs
###############################################################################

# mapResCode ()
#
sub mapResCode {
    local ($incode, $silent) = @_;
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
                         'GLU' => 'E'
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
#       &abort ("unknown residue '$incode'");
#       print STDERR ("unknown residue '$incode' (mapping to 'Z')\n");
#       $newcode = 'Z';
        if (! $silent) {
            print STDERR ("unknown residue '$incode' (mapping to 'X')\n");
        }
        $newcode = 'X';
    }
    elsif ($newcode eq 'X') {                  
        if (! $silent) {
            print STDERR ("strange residue '$incode' (seen code, mapping to 'X')
\n");
        }
    }

    return $newcode;
}


# coord_dist_sq ()
#
sub coord_dist_sq {
    my ($coords_a, $coords_b) = @_;
    my $x_delta = $coords_a->[0] - $coords_b->[0];
    my $y_delta = $coords_a->[1] - $coords_b->[1];
    my $z_delta = $coords_a->[2] - $coords_b->[2];
    return ($x_delta*$x_delta + $y_delta*$y_delta + $z_delta*$z_delta);
}


# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{usage: $0
\t -zonesfile     <zones_file>
\t -ssafile       <ssa_file>
\t -completefile  <complete_file>
\t -scoredfile    <scored_file>
\t[-outfile       <out_file>]
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "zonesfile=s",
		 "ssafile=s",
		 "completefile=s",
		 "scoredfile=s",
		 "outfile=s");

    # Check for legal invocation
    if (! defined $opts{zonesfile}    ||
	! defined $opts{ssafile} ||
	! defined $opts{completefile} ||
	! defined $opts{scoredfile}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{zonesfile});
    &checkExist ('f', $opts{ssafile});
    &checkExist ('f', $opts{completefile});
    &checkExist ('f', $opts{scoredfile});

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


# chip (chop for front of strings)
#
sub chip {
    my @flo = ();
    for ($i=0; $i <= $#_; ++$i) {
        $flo[$i] = substr ($_[$i], 0, 1);
        $_[$i] = substr ($_[$i], 1);                   # don't think this works
    }
    return $flo[0]  if ($#_ == 0);
    return @flo;
}


# chimp (chomp for front of strings)
#
sub chimp {
    my @flo = ();
    for ($i=0; $i <= $#_; ++$i) {
        $_[$i] =~ s/^(\s*)//;                          # don't think this works
        $flo[$i] = $1;
    }
    return $flo[0]  if ($#_ == 0);
    return @flo;
}


# cleanStr ()
#
sub cleanStr {
    my $str = shift;
    $str =~ s/[\x00-\x08\x0B-\x1F\x80-\xFF]//g;
    return $str;
}


# listMember ()
#
sub listMember {
    my ($item, @list) = @_;
    my $element;
    foreach $element (@list) {
        return $item  if ($item eq $element);
    }
    return undef;
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
    print "[$date][RUN][$0] $cmd\n" if ($debug && ! $silent);
    $ret = system ($cmd);
    #$ret = ($?>>8)-256;
    if ($ret != 0) {
	$date = `date +'%Y-%m-%d_%T'`;  chomp $date;
	print STDERR ("[$date][FAILURE:$ret][$0] $cmd\n");
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
	select (STDOUT);
    }
    print "[$date][LOG][$0] $msg\n";
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
            &alert ("dirnotfound: $path");
            exit -3;
	}
    }
    elsif ($type eq 'f') {
	if (! -f $path) {
            &alert ("filenotfound: $path");
            exit -3;
	}
	elsif (! -s $path) {
            &alert ("emptyfile: $path");
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
	&abort ("unable to open file $file for writing");
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
	    &abort ("unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("unable to open file $file for reading");
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
	    &abort ("unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("unable to open file $file for reading");
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
            &abort ("unable to open file $file for gzip -dc");
        }
    }
    elsif (! open (FILE, $file)) {
        &abort ("unable to open file $file for reading");
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
