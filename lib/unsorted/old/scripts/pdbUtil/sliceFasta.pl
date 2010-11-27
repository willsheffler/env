#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.1.1.1 $
##  $Date: 2003/09/05 01:47:28 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering

###############################################################################
# init
###############################################################################

if ($#ARGV < 2) {
    print STDERR "usage: $0 <fastafile> <startres[,startres2...]> <stopres[,stopres2...]>\n";
    exit -1;
}
$fastafile    = shift @ARGV;
$startres_opt = shift @ARGV;
$stopres_opt  = shift @ARGV;

# read fasta_in
@fasta_in = ();
foreach $line (&fileBufArray ($fastafile)) {
    next if ($line =~ /^\s*\>/);
    $line =~ s/\s+//g;
    push (@fasta_in, split (//, $line));
}

# get ranges
if ($startres_opt eq '_') {
    $startres_opt = 1;
}
if ($stopres_opt eq '_') {
    $stopres_opt = $#fasta_in + 1;
}
@startres = split (/,/, $startres_opt);
@stopres  = split (/,/, $stopres_opt);
&abort ("must define equal number of startres and stopres")  if ($#startres != $#stopres);

@ranges = ();
@mask = ();
for ($range_i=0; $range_i <= $#startres; ++$range_i) {
    &abort ("startres $startres[$range_i] is not digit")  if ($startres[$range_i] !~ /^\d+$/);
    &abort ("stopres $stopres[$range_i] is not digit")  if ($stopres[$range_i] !~ /^\d+$/);
    push (@ranges, "$startres[$range_i]-$stopres[$range_i]");
    &abort ("stopres $stopres[$range_i] must be >= startres $startres[$range_i]")  if ($stopres[$range_i] < $startres[$range_i]);
    for ($res_i=$startres[$range_i]; $res_i <= $stopres[$range_i]; ++$res_i) {
        &abort ("overlapping ranges (at least at residue $res_i)")  if ($mask[$res_i-1]);
        $mask[$res_i-1] = 1;
    }
}
$header .= " " . join (",", @ranges);

###############################################################################
# main
###############################################################################

# build fasta out
@fasta_out = ();
for ($res_i=0; $res_i <= $#fasta_in; ++$res_i) {
    next if (! $mask[$res_i]);
    push (@fasta_out, $fasta_in[$res_i]);
}

# output
$header  = $fastafile;
$header  =~ s!.*\/!!;
$header  =~ s/\.fasta//;
$header  =~ s/\.txt//;
$header .= " " . join (",", @ranges);
$buf = ">$header\n";
for ($i=0; $i <= $#fasta_out; ++$i) {
    $buf .= $fasta_out[$i];
    $buf .= "\n"  if (($i+1) % 50 == 0);
}
$buf .= "\n"  if ($buf !~ /\n$/);
print $buf;

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
    if (!defined $opts{file}
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
