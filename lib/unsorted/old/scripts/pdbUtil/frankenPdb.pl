#!/usr/bin/perl
##
##
## Copyright 2000, University of Washington
##   This document contains private and confidential information and
##   its disclosure does not constitute publication.  All rights are
##   reserved by University of Washington, except those specifically 
##   granted by license.
##
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.1.1.1 $
##  $Date: 2003/09/05 01:47:28 $
##  $Author: dylan $
##
##
###############################################################################

###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering

###############################################################################
# init
###############################################################################

# argv
local %opts = &getCommandLineOptions ();
$base_pdb = $opts{basepdb}; 
$outfile  = $opts{outfile};
@regions  = @{$opts{regions_list}};
@pdbs     = @{$opts{pdbs_list}};

###############################################################################
# main
###############################################################################

# fill frank with base_pdb
#
@buf = &fileBufArray ($base_pdb);
$atom_i = -1;
for ($line_i=0; $line_i <= $#buf; ++$line_i) {
    if ($buf[$line_i] =~ /^ATOM/ || $buf[$line_i] =~ /^HETATM/) {
	++$atom_i;
	$frank[$atom_i] = $buf[$line_i];
    }
}

# grab portions from other pdbs
#
for ($pdb_i=0; $pdb_i <= $#pdbs; ++$pdb_i) {
    @buf = &fileBufArray ($pdbs[$pdb_i]);

    # check size
    $atom_i = -1;
    for ($line_i=0; $line_i <= $#buf; ++$line_i) {
	if ($buf[$line_i] =~ /^ATOM/ || $buf[$line_i] =~ /^HETATM/) {
	    ++$atom_i;
	}
    }
    if ($atom_i != $#frank) {
	print STDERR "$0: pdbs must all have the same number of atoms and hetatms (base: ".($#frank+1)." vs. pdb".($pdb_i+1).": ".($atom_i+1).")\n";
	exit -2;
    }

    # replace portion of frank
    ($start_res,$stop_res) = split (/-/, $regions[$pdb_i]);
    $atom_i = -1;
    for ($line_i=0; $line_i <= $#buf; ++$line_i) {
	if ($buf[$line_i] =~ /^ATOM/ || $buf[$line_i] =~ /^HETATM/) {
	    ++$atom_i;
	    $res_n = substr ($buf[$line_i], 22, 4); 
	    if ($res_n >= $start_res && $res_n <= $stop_res) {
		$frank[$atom_i] = $buf[$line_i];
	    }
	}
    }
}

# write frank
#
if ($outfile) {
    open   (OUT, '>'.$outfile);
    select (OUT);
}
print join ("\n", @frank)."\n";
if ($outfile) {
    close  (OUT);
    select (STDOUT);
}

# done
#
exit 0;

###############################################################################
# subs
###############################################################################

# getCommandLineOptions()
#
#  desc: get the command line options
#
#  args: none
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    local $usage = qq{usage: $0
\t -basepdb  <base_pdb>
\t -regions  <startres1-stopres1,startres2-stopres2,...>
\t -pdbs     <pdb1,pdb2,...>
\t[-outfile  <outfile>]                                   (def: STDOUT)
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "basepdb=s", "regions=s", "pdbs=s", "outfile=s");
    
    # Check for legal invocation
    #
    if (! defined $opts{basepdb} ||
	! defined $opts{regions} ||
	! defined $opts{pdbs}
	) {
        print STDERR "$usage\n";
        exit -1;
    }
    
    $opts{regions_list} = +[];
    $opts{pdbs_list}    = +[];
    
    push (@{$opts{regions_list}}, split (/,/, $opts{regions}));
    push (@{$opts{pdbs_list}},    split (/,/, $opts{pdbs}));
    
    &checkExist ('f', $opts{basepdb});
    foreach $pdb (@{$opts{pdbs_list}}) {
	&checkExist ('f', $pdb);
    }
    
    return %opts;
}

###############################################################################
# util
###############################################################################

# runCmd()
#
sub runCmd {
    local $cmd = shift;
    local $retcode = 0;
#    print "RUN: '$cmd'\n\n";
    if (system ($cmd) != 0) {
#	&abort ("FAIL: $cmd");
	print STDERR ("FAIL: $cmd");
	$retcode = $?>>8;
    }
    return $retcode;
}

# checkExist()
#
sub checkExist {
    local ($type, $path) = @_;
    if ($type eq 'f' && ! -f $path) {
	&abort ("filenotfound $path");
	
    }
    if ($type eq 'd' && ! -d $path) { 
	&abort ("dirnotfound $path");
    }
    return 'true';
}

# checkExistAndCreate()
#
sub checkExistAndCreate {
    local ($type, $path) = @_;
    if ($type eq 'f' && ! -f $path) {
	print "creating $path...\n";
	open (FILE, '>'.$path);
	close (FILE);
    }
    if ($type eq 'd' && ! -d $path) { 
	print "creating $path...\n";
	$mode = 0777  if (! $mode);
	mkdir ($path, $mode);
    }
    return 'true';
}

# abort()
#
sub abort {
    local $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}

# fileBufString()
#
sub fileBufString {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

# fileBufArray()
#
sub fileBufArray {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
} 

###############################################################################
# end
###############################################################################
