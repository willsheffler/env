#!/usr/bin/perl
##
## Copyright 2001, University of Washington
##   This document contains private and confidential information and
##   its disclosure does not constitute publication.  All rights are
##   reserved by University of Washington, except those specifically
##   granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.6 $
##  $Date: 2005/06/21 01:41:23 $
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

$profit  = "$src_dir/shareware/ProFit/profit";
&abort ("unable to find $profit")  if (! -e $profit);

###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
my $pdb1         = $opts{pdb1};
my $pdb2         = $opts{pdb2};
my $pdb2fit_file = $opts{pdb2fit};
my $zonesfile    = $opts{zonesfile};

@zonesfile_buf = ($zonesfile) ? &fileBufArray ($zonesfile) : ();

###############################################################################
# main
###############################################################################

# read zones
foreach $line (@zonesfile_buf) {
    if ($line =~ /^zone\s*(\d+)\s*-\s*(\d+)\s*:\s*(\d+)\s*-\s*(\d+)/) {
	$zones_buf .= "zone $1-$2:$3-$4\n";
    }
}

# profit cmds
$profit_cmds  = qq{atoms CA\n};
$profit_cmds .= qq{$zones_buf}  if ($zones_buf);
$profit_cmds .= qq{fit\n};                   # get the rmsd
if ($pdb2fit_file) {
    $profit_cmds .= qq{write $pdb2fit_file\n};   # write out the rotated pdb
}

# run profit
#$rmsd = `echo '$profit_cmds' | $profit $pdb1 $pdb2 | fgrep RMS`;
$cmd = "echo '$profit_cmds' | $profit $pdb1 $pdb2";
$rmsd = `$cmd`;
$ret = $?;
$ret = $?>>8;
if ($ret != 0) {
    $ret -= 256  if ($ret >= 128);
    $date = `date +'%Y-%m-%d_%T'`;  $date =~ s/^\s+|\s+$//g;
    print STDERR ("[$date][FAILURE:$ret][$0] $cmd\n");
}
$rmsd =~ s/.*RMS\:\s*([\d\.]+)\s*/$1/s;    
$rmsd =~ s/Writing coordinates\.\.\.\s*$//;

# display rmsd
print "CA_RMSD: $rmsd\n";

# done
exit 0;

###############################################################################
# subs
###############################################################################


###############################################################################
# util
###############################################################################

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{usage: $0
\t -pdb1       <pdb1>
\t -pdb2       <pdb2>
\t[-zonesfile  <zonesfile>]
\t[-pdb2fit    <pdb2fit>]
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, "pdb1=s", "pdb2=s", "zonesfile=s", "pdb2fit=s");

    # Check for legal invocation
    if (!defined $opts{pdb1} ||
	!defined $opts{pdb2}

        ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{pdb1});
    &checkExist ('f', $opts{pdb2});
    &checkExist ('f', $opts{zonesfile}) if ($opts{zonesfile});

    return %opts;
}

# runCmd
#
sub runCmd {
    my $cmd = shift;
    &abort ("failure running cmd: $cmd") if (system ($cmd) != 0);
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
