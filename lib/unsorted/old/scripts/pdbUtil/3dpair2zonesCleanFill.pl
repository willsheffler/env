#!/usr/bin/perl
##
## Copyright 2003, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.1 $
##  $Date: 2004/04/27 20:21:15 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering
$debug = 1;

$min_seg_len = 3;

###############################################################################
# init
###############################################################################

# argv
#my %opts = &getCommandLineOptions ();
#my $file = $opts{file};

if ($#ARGV < 0) {
    print STDERR "usage: $0 <3dpair_file>\n";
    exit -1;
}
$file = shift @ARGV;

$file_buf = &bigFileBufArray ($file);

###############################################################################
# main
###############################################################################

# get alignment from file
#
foreach $line (@$file_buf) {
    $line =~ s/^\s+|\s+$//g;
    if (! @q_fasta) {
	($q_file, $q_seq) = split (/\s+/, $line);
	@q_fasta = split (//, $q_seq);
    }
    elsif (! @p_fasta) {
	($p_file, $p_seq) = split (/\s+/, $line);
	@p_fasta = split (//, $p_seq);
    }

    if ($line =~ /^[^\s\:]+\s+\d+\s+/) {
	@info = split (/\s+/, $line);
	$ranges    = $info[10];
	$alignment = $info[11];
	last;
    }
}
&abort ("couldn't find alignment line")  if (! $alignment);


# get ranges
#
$ranges =~ s/^\[|\]$//g;
($q_range, $p_range) = split (/\-/, $ranges);
($q_beg, $q_end) = split (/\:/, $q_range);
($p_beg, $p_end) = split (/\:/, $p_range);


# turn alignment into zones
#
@a = split (//, $alignment);
@zones = ();
@seqs  = ();
$in_seg = undef;
$q_i = $q_beg-1;
$p_i = $p_beg-1;
$q2p_map = ();
for ($i=0; $i <= $#a; ++$i) {
    if ($a[$i] =~ /[A-Z]/) {
	++$q_i;
	++$p_i;
    }
    elsif ($a[$i] =~ /[a-z]/) {
	++$p_i;
    }
    elsif ($a[$i] =~ /\-/) {
	++$q_i;
    }
    else {
	&abort ("unknown character '$a[$i]' in alignment at position ".($i+1));
    }

    if ($a[$i] =~ /[A-Z]/) {
	$q_stop = $q_i;
	$p_stop = $p_i;
	if (! $in_seg) {
	    $in_seg = 'true';
	    $q_start = $q_i;
	    $p_start = $p_i;
	}
	if ($i == $#a) {
	    if ($q_stop - $q_start + 1 >= $min_seg_len) {
		for ($map_i=0; $map_i <= $q_stop - $q_start; ++$map_i) {
		    $q2p_map[$q_start-1+$map_i] = $p_start-1+$map_i;
		}
	    }
	}
    }
    else {
	next if (! $in_seg);
	if ($q_stop - $q_start + 1 >= $min_seg_len) {
	    for ($map_i=0; $map_i <= $q_stop - $q_start; ++$map_i) {
		$q2p_map[$q_start-1+$map_i] = $p_start-1+$map_i;
	    }
	}
	$in_seg = undef;
    }
}
&abort ("query end not achieved: $q_i != $q_end")   if ($q_i != $q_end);
&abort ("parent end not achieved: $p_i != $p_end")  if ($p_i != $p_end);


# fill in gaps of equal length
#
$in_seq = undef;
$last_def_qi = -1;
$last_def_pj = -1;
for ($qi=0; $qi <= $#q_fasta; ++$qi) {
    if (! defined $q2p_map[$qi]) {
	$next_def_qi = undef;
	for ($qi2=$qi+1; $qi2 <= $#q_fasta+1; ++$qi2) {
	    if ($qi2 == $#q_fasta+1) {
		$next_def_qi = $#q_fasta+1;
		$next_def_pj = $#p_fasta+1;
		last;
	    } elsif (defined $q2p_map[$qi2]) {
		$next_def_qi = $qi2;
		$next_def_pj = $q2p_map[$qi2];
		last;
	    }
	}
	if ($next_def_qi - $last_def_qi == $next_def_pj - $last_def_pj) {
	    for ($map_i=1; $map_i < $next_def_qi - $last_def_qi; ++$map_i) {
		$q2p_map[$last_def_qi+$map_i] = $last_def_pj+$map_i;
	    }
	}
    }
    else {
	$last_def_qi = $qi;
	$last_def_pj = $q2p_map[$qi];
    }
}
# do short termini
for ($qi=0; $qi <= $#q_fasta; ++$qi) {
    if (defined $q2p_map[$qi]) {    
	if ($qi > 0) {
	    $pj = $q2p_map[$qi];
	    if ($qi < 15 || $pj < 15) {
		for ($qi2=$qi-1; $qi2 >= 0; --$qi2) {
		    --$pj;
		    last if ($pj < 0);
		    $q2p_map[$qi2] = $pj;
		}
	    }
	}
	last;
    }
}
for ($qi=$#q_fasta; $qi >= 0; --$qi) {
    if (defined $q2p_map[$qi]) {    
	if ($qi < $#q_fasta) {
	    $pj = $q2p_map[$qi];
	    if (($q_fasta - $qi < 15) || ($p_fasta - $pj < 15)) {
		for ($qi2=$qi+1; $qi2 <= $#q_fasta; ++$qi2) {
		    ++$pj;
		    last if ($pj > $#p_fasta);
		    $q2p_map[$qi2] = $pj;
		}
	    }
	}
	last;
    }
}



# build zones
#
@zones = ();
for ($qi=0, $start_qi=undef, $last_pj=-1000; 
     $qi <= $#q2p_map+1; 
     ++$qi) {
    $pj = $q2p_map[$qi];
    if (! defined $pj || (defined $pj && $pj != $last_pj+1)) {
        if (defined $start_qi) {
            ++$start_qi;
            ++$start_pj;
            ++$last_qi;
            ++$last_pj;
	    push (@zones, sprintf (qq{zone %4d-%-4d:%4d-%-4d}, $start_qi, $last_qi, $start_pj, $last_pj));
        }
        if (defined $pj) {
            $start_qi = $qi;
            $start_pj = $pj;
        } else {
            $start_qi = undef;
            $start_pj = undef;
        }
    }
    $last_qi = $qi;
    $last_pj = $pj  if (defined $pj);
}



# output
#
$out = join ("\n", @zones)."\n";
#for ($i=0; $i <= $#zones; ++$i) {
##    $out .= $zones[$i] ."\t". $segs[$i] ."\n";
#    $out .= $zones[$i] ."\n";
#}
print $out;


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
    my $usage = qq{usage: $0 -file <file>\n};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, "file=s");

    # Check for legal invocation
    if (! defined $opts{file}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{file});

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

# abort()
#
sub abort {
    my $msg = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print STDERR "[$date]:$0:ABORT: $msg\n";
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
