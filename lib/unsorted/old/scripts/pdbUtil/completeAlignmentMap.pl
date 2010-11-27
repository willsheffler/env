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
##  $Date: 2005/04/20 20:25:49 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering
$debug = 1;

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

$ksync   = "$src_dir/ksync/bin/ksync"; 

###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
my $fasta_file     = $opts{fastafile};
my $alignment_file = $opts{alignmentfile};
my $out_file       = $opts{outfile};

###############################################################################
# main
###############################################################################

# ksync align incomplete sequence with complete sequence in fasta file
#
$incomplete_fasta = '';
foreach $line (&fileBufArray ("$alignment_file")) {
    next if ($line =~ /^>/);
    $line =~ s/[^a-zA-Z]//g;
    $incomplete_fasta .= $line;
}
$incomplete_fasta = uc $incomplete_fasta;
$incomplete_fasta_file = "$tmp_dir/incomplete.fasta";
open (INCFASTA, '>'.$incomplete_fasta_file);
print INCFASTA $incomplete_fasta."\n";
close (INCFASTA);
$incomplete2completeMap_zones_file = "$tmp_dir/incomplete2completeMap.zones";
&runCmd (qq{$ksync -r l -v_i .001 -v_e .0001 -m f,_,_,_,_,_,_,_,_,_,_,_,_ -a $incomplete_fasta_file -A $fasta_file -d 2 -f_ksync $incomplete2completeMap_zones_file}, undef, 'SILENT');
unlink $incomplete_fasta_file;


# read in incomplete2completeMap
#
@incomplete2completeMap = ();
foreach $line (&fileBufArray ($incomplete2completeMap_zones_file)) {
    next if ($line !~ s/^zone\s*//i);
    ($start_1, $stop_1, $start_2, $stop_2) = split (/\D+/, $line);
    &abort ("unequal zones $line")  if ($stop_1-$start_1 != $stop_2-$start_2);
    for ($i=0; $i <= $stop_1-$start_1; ++$i) {
	$incomplete2completeMap[$start_1+$i-1] = $start_2+$i-1;
    }
}
unlink $incomplete2completeMap_zones_file;


# build realMap
#
@realMap = ();
$query_i = -1;
$incomplete_i = -1;
foreach $line (&fileBufArray ("$alignment_file")) {
    next if ($line =~ /^>/);
    $line =~ s/[^a-zA-Z\.\-]//g;
    foreach $c (split (//, $line)) {
	++$query_i;
	if ($c =~ /[a-zA-Z]/) {
	    ++$incomplete_i;
	    $realMap[$query_i] = $incomplete2completeMap[$incomplete_i];
	}
    }
}
# DEBUG
#for ($i=0; $i <= $#realMap; ++$i) {
#    print "$i -> $realMap[$i]\n";
#}


# build output zones
#
@out = ();
for ($qi=0, $start_qi=undef, $last_pj=undef; 
     $qi <= $#realMap+1; 
     ++$qi) {
    $pj = $realMap[$qi];
    if (! defined $pj || 
	(defined $pj && ! defined $last_pj) ||
	(defined $pj && defined $last_pj && $pj != $last_pj+1) 
	) {
        if (defined $start_qi) {
            ++$start_qi;
            ++$start_pj;
            ++$last_qi;
            ++$last_pj;
	    push (@out, sprintf (qq{zone %4d-%-4d:%4d-%-4d}, $start_qi, $last_qi, $start_pj, $last_pj));
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
    $last_pj = $pj;
}


# output
if ($out_file) {
    open (OUT, '>'.$out_file);
    select OUT;
}
print join ("\n", @out)."\n";
if ($out_file) {
    close (OUT);
    select STDOUT;
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
    my $usage = qq{usage: $0
\t -fastafile      <fasta_file>
\t -alignmentfile  <alignment_file>
\t[-outfile        <out_file>]
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, "fastafile=s", "alignmentfile=s", "outfile=s");

    # Check for legal invocation
    if (! defined $opts{fastafile} ||
	! defined $opts{alignmentfile}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{fastafile});
    &checkExist ('f', $opts{alignmentfile});

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
