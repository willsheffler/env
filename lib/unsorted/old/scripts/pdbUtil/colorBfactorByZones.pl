#!/usr/bin/perl
##
## Copyright 2004, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.1 $
##  $Date: 2004/07/19 00:10:34 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$| = 1;
$debug = 1;

###############################################################################
# init
###############################################################################

if ($#ARGV < 2) {
    print STDERR "usage: $0 <pdb> <outfile> <zonefile> [q/p(def:q)]\n";
    exit -1;
}
$rosetta_pdb = shift @ARGV;
$out_file    = shift @ARGV;
$zonefile    = shift @ARGV;
$q_or_p      = shift @ARGV;
$q_or_p = 'q'  if (! $q_or_p);

&checkExist ('f', $rosetta_pdb);
&checkExist ('f', $zonefile);

@pdb_buf   = &fileBufArray ($rosetta_pdb);
@zones_buf = &fileBufArray ($zonefile);

###############################################################################
# main
###############################################################################

# read zones
#
foreach $line (@zones_buf) {
    next if ($line !~ /^\s*zone/i);
    $line =~ s/^\D+//g;
    @zones = split (/\D+/, $line);
    if ($q_or_p =~ /q/i) {
	for ($res_i=$zones[0]; $res_i <= $zones[1]; ++$res_i) {
	    $template[$res_i] = 'true';
	}
    } else {
	for ($res_i=$zones[2]; $res_i <= $zones[3]; ++$res_i) {
	    $template[$res_i] = 'true';
	}
    }
}


# body
#
@out = ();
for ($i=0; $i <= $#pdb_buf; ++$i) {
    if ($pdb_buf[$i] =~ /^ATOM|^HETATM/) {

	# fix any missing occupancy and temp factor field
	$linelen  = length ($pdb_buf[$i]);
	$pdb_buf[$i] .= ' 'x(79-$linelen);
	$prev_occup = substr ($pdb_buf[$i], 54, 6);
	$prev_bfact = substr ($pdb_buf[$i], 60, 6);
	if ($prev_occup eq ' 'x6) {
	    substr ($pdb_buf[$i], 54, 6) = sprintf ("%6.2f", 1.00);
	}
	if ($prev_bfact eq ' 'x6) {
	    substr ($pdb_buf[$i], 60, 6) = sprintf ("%6.2f", 0.00);
	}

	# get res number residue numbering if necessary
	$res_n = substr ($pdb_buf[$i], 22, 4);

	# adjust temperature factor to include coloring
	if ($pdb_buf[$i] =~ /^HETATM/) {
	    substr ($pdb_buf[$i], 60, 6) = sprintf ("%6.2f", 4.00);
	}
	elsif ($template[$res_n]) {
	    substr ($pdb_buf[$i], 60, 6) = sprintf ("%6.2f", 1.00);
	}
	else {
	    substr ($pdb_buf[$i], 60, 6) = sprintf ("%6.2f", 10.00);
	}
    }
    push (@out, $pdb_buf[$i]);
}


# write it
#
open (OUT, '>'.$out_file);
print OUT join ("\n", @out) . "\n";
close(OUT);


# done
exit 0;

###############################################################################
# util
###############################################################################

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
    $ret = ($?>>8);
    if ($ret != 0) {
	$ret -= 256  if ($ret >= 128);
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
###############################################################################
