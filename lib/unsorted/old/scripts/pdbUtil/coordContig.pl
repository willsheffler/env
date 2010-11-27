#!/usr/bin/perl -w
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
##  $Revision: 1.2 $
##  $Date: 2004/03/31 23:11:49 $
##  $Author: dylan $
##
##
###############################################################################


###############################################################################
# init
###############################################################################

$| = 1;                                                   # don't buffer stdout

local %opts = &getCommandLineOptions ();
@pdbfile = &fileBufArray ($opts{pdbfile});
$chainOI = $opts{chain};

###############################################################################
# main
###############################################################################

foreach $line (@pdbfile) {
    if ($line =~ /^ATOM/) {
	$chain = substr ($line, 21, 1);
	if ($chain eq $chainOI) {
	    $resSeq = substr ($line, 22, 4);
	    $resName = substr ($line, 17, 3);
	    $occ     = substr ($line, 54, 6);
	    if ($resSeq >= 0 && $occ >= 0) {
		$residues[$resSeq] = $resName;
	    }
	}
    }
}
if (! @residues) {
    print "no such chain '$chainOI' found\n";
    exit -1;
}


$in_seg = undef;
for ($i=1; $i <= $#residues+1; ++$i) {
    if (! $residues[$i] && $in_seg) {
        $close = $i-1;
        print "contig: $residues[$open]$open-$residues[$close]$close\n";
        $in_seg = undef;
    }
    elsif ($residues[$i] && ! $in_seg) {
        $in_seg = 'true';
        $open = $i;
    }
}
 
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
    local $usage = qq{usage: $0 -pdbfile <pdbfile> -chain <chain>};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "pdbfile=s", "chain=s");

    # Check for legal invocation
    #
    if (! defined $opts{'pdbfile'}) {
        print STDERR "$usage\n";
        exit -1;
    }

    # defaults
    #
    if ($opts{chain}) {
	$opts{chain} = &chip ($opts{chain});
	$opts{chain} = ($opts{chain} eq '_') ? ' ' : uc $opts{chain};
    } else {
	$opts{chain} = ' ';
    }

    return %opts;
}

###############################################################################
# util
###############################################################################

sub chip {
    my @flo = ();
    for ($i=0; $i <= $#_; ++$i) {
	$flo[$i] = substr ($_[$i], 0, 1);
	$_[$i] = substr ($_[$i], 1);                   # don't think this works
    }
    return $flo[0]  if ($#_ == 0);
    return @flo;
}

sub chimp {
    my @flo = ();
    for ($i=0; $i <= $#_; ++$i) {
	$_[$i] =~ s/^(\s*)//;	                       # don't think this works
	$flo[$i] = $1;
    }
    return $flo[0]  if ($#_ == 0);
    return @flo;
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
