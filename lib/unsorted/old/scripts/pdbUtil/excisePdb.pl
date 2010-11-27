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
    print STDERR "usage: $0 <fullpdb> <chain> <startres> <stopres> <outfile>\n";
    exit -1;
}
$pdb      = shift @ARGV;
$chainOI  = shift @ARGV;
$startres = shift @ARGV;
$stopres  = shift @ARGV;
$outfile  = shift @ARGV;

if (! -f $pdb) {
    print STDERR "$0: filenotfound: $pdb\n";
    exit -1;
}
@pdb_buf = &fileBufArray ($pdb);

$chainOI  = ($chainOI eq '_')  ? ' '      : uc $chainOI;
$startres = ($startres eq '_') ? -1000000 : $startres;
$stopres  = ($stopres eq '_')  ?  1000000 : $stopres;

if ($startres =~ s/([a-z])$//i) {
    $startres_iCode = $1;
}
if ($stopres =~ s/([a-z])$//i) {
    $stopres_iCode = $1;
}

# body
$out_buf = '';
$res_cnt = 0;
for ($i=0; $i <= $#pdb_buf; ++$i) {
    last if ($chainOI_found && ($pdb_buf[$i] =~ /^MODEL/ || $pdb_buf[$i] =~ /^ENDMDL/)); 
    if ($pdb_buf[$i] =~ /^ATOM|^HETATM/) {
	$header_done = 'TRUE';
	$altLoc = substr ($pdb_buf[$i], 16, 1);
	$chain  = substr ($pdb_buf[$i], 21, 1);
	$resseq = substr ($pdb_buf[$i], 22, 4);
	$iCode  = substr ($pdb_buf[$i], 26, 1);
	$iCode  =~ s/\s+//; 
	next if ($chain ne $chainOI);
	next if ($altLoc ne ' ' && $altLoc ne 'A' && $altLoc ne '1');
	#next if ($resseq < $startres || $resseq > $stopres);
	#next if ($resseq eq $startres && $iCode ne $startres_iCode);
	#next if ($resseq ne $stopres && $iCode ne $stopres_iCode);

	next if ($resseq >= $startres && $resseq <= $stopres);
	$chainOI_found = 'true';

	if ($resseq ne $last_resseq) {
	    ++$res_cnt;
	    $last_resseq = $resseq;
	}
	last if ($stopped && ($resseq != $stopres || $iCode ne $stopres_iCode));
	$out_buf .= $pdb_buf[$i]."\n";
    }
    $out_buf .= $pdb_buf[$i]."\n"  if (! $header_done);
}

if ($chainOI_found) {
    open  (OUT, '>'.$outfile);
    print  OUT $out_buf;
    close (OUT);
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
