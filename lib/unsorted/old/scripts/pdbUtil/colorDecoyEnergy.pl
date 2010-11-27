#!/usr/bin/perl
##
## Copyright 2003, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.3 $
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

###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
my $decoy_file = $opts{decoyfile};
my $field      = $opts{field};
my $avg        = $opts{avg};
my $out_file   = $opts{outfile};

#if ($#ARGV < 0) {
#    print STDERR "usage: $0 <file>\n";
#    exit -1;
#}
#$file = shift @ARGV;

$decoy_buf = &bigFileBufArray ($decoy_file);

###############################################################################
# main
###############################################################################

# ROCK ON
$started = undef;
$started_avg = undef;
$field_i = undef;
@energies = ();
foreach $line (@$decoy_buf) {
    if ($avg && $line =~ /^\s*energies-average/) {
	$started_avg = 'true';
	next;
    }
    next if ($avg && ! $started_avg);

    if ($line =~ /^\s*res/) {
	$line =~ s/^\s+|\s+$//g;
	@fields = split (/\s+/, $line);
	for ($header_i=0; $header_i <= $#fields; ++$header_i) {
	    if ($fields[$header_i] eq $field) {
		$field_i = $header_i;
		last;
	    }
	}
	&abort ("can't find field $field")  if (! defined $field_i);
	$started = 'true';
	next;
    }

    next if (! $started);
    last if ($line !~ /^\s*\d+\s*/);

    $line =~ s/^\s+|\s+$//g;
    @scores = split (/\s+/, $line);
    push (@energies, $scores[$field_i]);
}


# get spread
#
$min_energy =  1000000;
#$max_energy = -1000000;
foreach $energy (@energies) {
    $min_energy = $energy  if ($energy < $min_energy);
#    $max_energy = $energy  if ($energy > $max_energy);
}
#$delta = $max_energy - $min_energy;


# determine Bvalues
#
@Bvals = ();
for ($res_i=0; $res_i <= $#energies; ++$res_i) {
    $Bvals[$res_i] = 10 * ($energies[$res_i] - $min_energy) + 1;
}

# adjust Bfactors
$res_i = -1;
$last_resseq = undef;
for ($line_i=0; $line_i <= $#{$decoy_buf}; ++$line_i) {
    if ($decoy_buf->[$line_i] =~ /^ATOM  /) {
	$resseq = substr ($decoy_buf->[$line_i], 22, 4);
	if (! defined $last_resseq || $last_resseq != $resseq) {
	    ++$res_i;
	    $last_resseq = $resseq;
	}

	substr ($decoy_buf->[$line_i], 60, 6) = sprintf ("%6.2f", $Bvals[$res_i]);
    }
}

if ($out_file) {
    open (OUT, '>'.$out_file);
    select (OUT);
}
print join ("\n", @{$decoy_buf}) ."\n"; 
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
    my $usage = qq{usage: $0 -decoyfile <decoy_file> [-field <field(e.g.Eatr/Erep/Eres,def:Eres)>] [-avg] [-outfile <out_file>]\n};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, "decoyfile=s", "field=s", "avg", "outfile=s");

    # Check for legal invocation
    if (! defined $opts{decoyfile}
	
        ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{decoyfile});

    $opts{field} = 'Eres'  if (! defined $opts{field});

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
