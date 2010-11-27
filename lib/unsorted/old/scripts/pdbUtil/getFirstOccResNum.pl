#!/usr/bin/perl -w
##
##
## Copyright 2001, University of Washington
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
package PDButil;
###############################################################################

###############################################################################
# init
###############################################################################

$| = 1;                                                   # don't buffer stdout

local %opts = &getCommandLineOptions ();
local $pdbfile   = $opts{pdbfile};
local $chainOI   = $opts{chain};
local $searchres = $opts{searchres}; 
local $direction = $opts{direction};

###############################################################################
# main
###############################################################################

# read
#
@buf = &fileBufArray ($pdbfile);

$res_num = 0;
for ($i=0; $i <= $#buf; ++$i) {
    last if ($chainOI_found && ($buf[$i] =~ /^TER/ || $buf[$i] =~ /^MODEL/ || $buf[$i] =~ /^ENDMDL/));
    last if ($chainOI_found && $buf[$i] =~ /^T/);
    next if ($buf[$i] !~ /^ATOM/);# && $buf[$i] !~ /^HETATM/);
    $chain = substr ($buf[$i], 21, 1);
    next if ($chain ne $chainOI);
    $chainOI_found = 'TRUE';
    $altLoc = substr ($buf[$i], 16, 1);
    next if ($altLoc ne ' ' && $altLoc ne 'A');
    $atomtype = substr ($buf[$i], 12, 4);
    next if ($atomtype ne ' CA ');

    ++$res_num;
    $occ = substr ($buf[$i], 54, 6);
    if ($direction eq 'first') {
	if ($res_num >= $searchres) {
	    if ($occ > 0) {
		print "$res_num\n";
		exit 0;
	    }
	}
    }
    elsif ($direction eq 'last') {
	if ($occ > 0) {
	    $last_occ = $res_num;
	}
	if ($searchres) {
	    if ($res_num >= $searchres) {
		print "$last_occ\n";
		exit 0;
	    }
	}
    }
    else {
	&abort ("unknown direction $direction");
    }
}
if (! $chainOI_found) {
    &abort ("no such chain found '$chainOI' in pdb file $pdbfile");
}

# in case the TER wasn't found
if ($direction eq 'last') {
    print "$last_occ\n";
    exit 0;
}

# should never be here
print STDERR "FAILURE IN $0 LOGIC\n";
exit -1;

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
\t -pdbfile   <pdbfile>
\t[-chain     <chain>]       (def: _)
\t[-searchres <searchres>]   (def: off)
\t[-direction <first/last>]  (def: first)
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "pdbfile=s", "chain=s", "searchres=i", "direction=s");


    # Check for legal invocation
    #
    if (! defined $opts{pdbfile}
	) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExistence ('f', $opts{pdbfile});	


    # Defaults
    #
    $opts{chain}     = ' '      if (! defined $opts{chain} || 
				    $opts{chain} eq '0' || 
				    $opts{chain} eq '_');
    $opts{chain}     = uc $opts{chain};
    $opts{direction} = 'first'  if (! defined $opts{direction});
    $opts{searchres} = 1        if (! defined $opts{searchres} && 
				    $opts{direction} eq 'first');
    $opts{searchres} = undef    if (! defined $opts{searchres} && 
				    $opts{direction} eq 'last');

 
    return %opts;
}

###############################################################################
# util
###############################################################################

sub maxInt {
    local ($v1, $v2) = @_;
    return ($v1 > $v2) ? $v1 : $v2;
}

sub tidyDecimals {
    my ($num, $decimal_places) = @_;
    if ($num !~ /\./) {
	$num .= '.' . '0' x $decimal_places;
	$num =~ s/^0+//;
    }
    else {
	if ($num =~ s/(.*\.\d{$decimal_places})(\d).*$/$1/) {
	    my $nextbit = $2;
	    if ($nextbit >= 5) {
		my $flip = '0.' . '0' x ($decimal_places - 1) . '1'; 
		$num += $flip;
	    }
        }
	$num =~ s/^0//;
	my $extra_places = ($decimal_places + 1) - length $num;
	$num .= '0' x $extra_places  if ($extra_places > 0);
    }

    return $num;
}

sub distsq {
    local @dims = @_;
    local $v = 0;
    foreach $dim (@dims) {
	$v += $dim*$dim;
    }
    return $v;
}

sub logMsg {
    local ($msg, $logfile) = @_;

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    print $msg, "\n";
    if ($logfile) {
        close (LOGFILE);
        select (STDOUT);
    }
    return 'true';
}

sub checkExistence {
    local ($type, $path) = @_;
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
    }
}

sub abort {
    local $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}

sub writeBufToFile {
    ($file, $bufptr) = @_;
    if (! open (FILE, '>'.$file)) {
	&abort ("$0: unable to open file $file for writing");
    }
    print FILE join ("\n", @{$bufptr}), "\n";
    close (FILE);
    return;
}

sub fileBufString {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

sub fileBufArray {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
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
1; # package end
# end
###############################################################################
