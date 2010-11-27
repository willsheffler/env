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
                   

# init
if ($#ARGV < 1) {
    print STDERR "usage: $0 <pdb> <Bvalue> [<atom/hetero/both(def:both)>]\n";
    exit -1;
}
$pdb            = shift @ARGV;
$Bvalue         = shift @ARGV;
$atom_or_hetero = shift @ARGV;
$atom_or_hetero = (! $atom_or_hetero) ? 'both' : $atom_or_hetero;

# body
if (! -f $pdb) {
    print STDERR "$0: filenotfound: $pdb\n";
    exit -1;
}
@pdb_buf = &fileBufArray ($pdb);
for ($i=0; $i <= $#pdb_buf; ++$i) {
    if ($pdb_buf[$i] =~ /^ATOM/ && ($atom_or_hetero =~ /^a/i ||
				    $atom_or_hetero =~ /^b/i)) {

	# fix any missing occupancy and temp factor field
	$linelen  = length ($pdb_buf[$i]);
	$pdb_buf[$i] .= ' 'x(79-$linelen);
	$prev_occup = substr ($pdb_buf[$i], 54, 6);
	$prev_bfact = substr ($pdb_buf[$i], 60, 6);
	if ($prev_occup eq ' 'x6) {
	    substr ($pdb_buf[$i], 54, 6) = sprintf ("%6.2f", 1);
	}
	substr ($pdb_buf[$i], 60, 6) = sprintf ("%6.2f", $Bvalue);
    }
    elsif ($pdb_buf[$i] =~ /^HETATM/ && ($atom_or_hetero =~ /^h/i ||
					 $atom_or_hetero =~ /^b/i)) {

	# fix any missing occupancy and temp factor field
	$linelen  = length ($pdb_buf[$i]);
	$pdb_buf[$i] .= ' 'x(79-$linelen);
	$prev_occup = substr ($pdb_buf[$i], 54, 6);
	$prev_bfact = substr ($pdb_buf[$i], 60, 6);
	if ($prev_occup eq ' 'x6) {
	    substr ($pdb_buf[$i], 54, 6) = sprintf ("%6.2f", 1);
	}
	substr ($pdb_buf[$i], 60, 6) = sprintf ("%6.2f", $Bvalue);
    }
    print $pdb_buf[$i] , "\n";
}

exit 0;


###############################################################################
# util
###############################################################################

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
