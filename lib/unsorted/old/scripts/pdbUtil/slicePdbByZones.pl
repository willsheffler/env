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
#
if ($#ARGV < 1) {
    print STDERR "usage: $0 <zonesfile> <query/parent> <fullpdb> <chain> <outfile> [renum_flag] [sequential_flag]\n";
    exit -1;
}
$zonesfile       = shift @ARGV;
$q_or_p          = shift @ARGV;
$pdb             = shift @ARGV;
$chainOI         = shift @ARGV;
$outfile         = shift @ARGV;
$renum_flag      = shift @ARGV  if (@ARGV);
$sequential_flag = shift @ARGV  if (@ARGV);

if (! -f $pdb) {
    print STDERR "$0: filenotfound: $pdb\n";
    exit -1;
}
@pdb_buf = &fileBufArray ($pdb);

$chainOI  = ($chainOI eq '_')  ? ' '      : uc $chainOI;


@startres = ();
@startres = ();
foreach $line (&fileBufArray ($zonesfile)) {
    next if ($line !~ /^\s*zone\s*(\d+)\s*-\s*(\d+)\s*:\s*(\d+)\s*-\s*(\d+)/i);
    $q_start = $1;
    $q_stop  = $2;
    $p_start = $3;
    $p_stop  = $4;
    
    if ($q_or_p =~ /^p/i) {
	push (@startres, $p_start);
	push (@stopres, $p_stop);
    }
    else {
	push (@startres, $q_start);
	push (@stopres, $q_stop);
    }
}


@mask = ();
for ($range_i=0; $range_i <= $#startres; ++$range_i) {
    if ($startres[$range_i] =~ s/([a-z])$//i) {
	$startres_iCode[$range_i] = $1;
    }
    if ($stopres[$range_i] =~ s/([a-z])$//i) {
	$stopres_iCode[$range_i] = $1;
    }
    &abort ("stopres $stopres[$range_i] must be >= startres $startres[$range_i]")  if ($stopres[$range_i] < $startres[$range_i]);
    for ($res_i=$startres[$range_i]; $res_i <= $stopres[$range_i]; ++$res_i) {
	&abort ("overlapping ranges (at least at residue $res_i)")  if ($mask[$res_i]);
	$mask[$res_i] = 1;
    }
}

# body
$out_buf = '';
$res_cnt  = 0;
$atom_cnt = 0;
$something_found = undef;
for ($range_i=0; $range_i <= $#startres; ++$range_i) {
    $chainOI_found = undef;
    $header_done   = undef;
    $started       = undef;
    $stopped       = undef;
    $last_resseq   = -100000;
    foreach $line (@pdb_buf) {
	last if ($chainOI_found && ($line =~ /^TER/ || $line =~ /^MODEL/ || $line =~ /^ENDMDL/)); 
	if ($line =~ /^ATOM|^HETATM/) {
	    $header_done = 'TRUE';
	    $altLoc = substr ($line, 16, 1);
	    $chain  = substr ($line, 21, 1);
	    $resseq = substr ($line, 22, 4);
	    $iCode  = substr ($line, 26, 1);
	    $iCode  =~ s/\s+//; 
	    next if ($chain ne $chainOI);
	    next if ($altLoc ne ' ' && $altLoc ne 'A' && $altLoc ne '1');
	    #next if ($resseq < $startres || $resseq > $stopres);
	    #next if ($resseq eq $startres && $iCode ne $startres_iCode);
	    #next if ($resseq ne $stopres && $iCode ne $stopres_iCode);
	    if ($startres[$range_i] == -1000000 || $resseq > $startres[$range_i] || ($resseq == $startres[$range_i] && $iCode eq $startres_iCode[$range_i])) {
		$started = 'TRUE';
	    }
	    next if (! $started);
	    if ($resseq > $stopres[$range_i]) {
		last;
	    }
	    if ($resseq == $stopres[$range_i] && $iCode eq $stopres_iCode[$range_i]) {
		$stopped = 'TRUE';
	    }
	    $something_found = 'TRUE';
	    $chainOI_found = 'TRUE';
	    if ($resseq ne $last_resseq) {
		++$res_cnt;
		$last_resseq = $resseq;
	    }
	    ++$atom_cnt;
	    if ($renum_flag) {
		substr ($line,  6, 5) = sprintf ("%5d", $atom_cnt);
		if ($sequential_flag) {
		    substr ($line, 22, 4) = sprintf ("%4d", $res_cnt);
		} else {
		    substr ($line, 22, 4) = sprintf ("%4d", $resseq-$startres[0]+1);
		}
	    }
	    last if ($stopped && ($resseq != $stopres[$range_i] || $iCode ne $stopres_iCode[$range_i]));
	    $out_buf .= $line."\n";
	}
    }
}

if ($something_found) {
    if ($outfile) {
	open (OUT, '>'.$outfile);
	select (OUT);
    }
    print $out_buf;
    if ($outfile) {
	close (OUT);
	select (STDOUT);
    }
}
else {
    &abort ("no density found within range $startres:$stopres for chain $chainOI in file $pdb");
}

exit 0;


###############################################################################
# util
###############################################################################

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
# end
###############################################################################
