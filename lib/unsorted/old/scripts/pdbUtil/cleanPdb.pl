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
##  $Revision: 1.1.1.1 $
##  $Date: 2003/09/05 01:47:28 $
##  $Author: dylan $
##
##
###############################################################################

###############################################################################
package PDButil;
###############################################################################

$| = 1;                                                   # don't buffer stdout

local %opts = &getCommandLineOptions ();
local $pdbfile = $opts{pdbfile};
local $outfile = $opts{outfile};

$pdbID = $pdbfile;
$pdbID =~ s!^.*/?p?d?b?(\w\w\w\w)\.[pe][dn][bt]\.?g?z?Z?$!$1!;
$pdbID = lc $pdbID;

###############################################################################
# main
###############################################################################

# read
#
#@pdbbuf = &fileBufArray ($pdbfile);
#foreach $line (@pdbbuf) {
#    if ($line =~ /^ATOM/) {
#	$occ = substr ($line, 54, 6);
#	if ($occ > 0.0) {
#	    push (@out, $line);
#	}
#    } elsif ($line !~ /^HETATM/) {
#	push (@out, $line);
#    }
#}

@pdbbuf = &fileBufArray ($pdbfile);
$last_chain  = '';
$last_resseq = undef;
@aa_buf = ();
$bb_defined = undef;
for ($i=0; $i <= $#pdbbuf; ++$i) {
    $line = $pdbbuf[$i];
    if ($line =~ /^ATOM/) {

	$occ = substr ($line, 54, 6);
	next if ($occ <= 0.0);

	$chain = substr ($line, 21, 1);
	if ($chain ne $last_chain) {
	    $last_chain = $chain;
	    $last_resseq = undef;
	    @aa_buf = ();
	    $N_occ      = undef;
	    $CA_occ     = undef;
	    $C_occ      = undef;
	    $O_occ      = undef;
	    $bb_defined = undef;
	}

	$resseq = substr ($line, 22, 4);
	$last_resseq = $resseq  if (! defined $last_resseq);
	if ($resseq != $last_resseq || $i == $#pdbbuf) {
	    if ($i == $#pdbbuf) {
		substr ($line, 21, 1) = ' ';                     # remove chain
		push (@aa_buf, $line);
		$atomname = substr ($line, 12, 4);
		if ($atomname eq ' N  ') {
		    $N_occ  = 1;
		} elsif ($atomname eq ' CA ') {
		    $CA_occ = 1;
		} elsif ($atomname eq ' C  ') {
		    $C_occ  = 1;
		} elsif ($atomname eq ' O  ') {
		    $O_occ  = 1;
		}
		$bb_defined = 'true'  if ($N_occ && $CA_occ && $C_occ && $O_occ);
	    }
	    if ($bb_defined) {
		$bb_defined_somewhere = 'true';
		push (@out, @aa_buf);
	    } else {
		print "DISCARDING chain:'$chain' resseq:'$resseq'\n";
	    }
	    $last_resseq = $resseq;
	    @aa_buf = ();
	    $N_occ      = undef;
	    $CA_occ     = undef;
	    $C_occ      = undef;
	    $O_occ      = undef;
	    $bb_defined = undef;
	}
	substr ($line, 21, 1) = ' ';                             # remove chain
	push (@aa_buf, $line);
	$atomname = substr ($line, 12, 4);
	if ($atomname eq ' N  ') {
	    $N_occ  = 1;
	} elsif ($atomname eq ' CA ') {
	    $CA_occ = 1;
	} elsif ($atomname eq ' C  ') {
	    $C_occ  = 1;
	} elsif ($atomname eq ' O  ') {
	    $O_occ  = 1;
	}
	$bb_defined = 'true'  if ($N_occ && $CA_occ && $C_occ && $O_occ);
    }
    else {
	if (@aa_buf) {
	    if ($bb_defined) {
		push (@out, @aa_buf);
		$bb_defined_somewhere = 'true';
	    } else {
		print "DISCARDING chain:'$chain' resseq:'$resseq'\n";
	    }
	    $last_chain  = undef;
	    $last_resseq = undef;
	    @aa_buf = ();
	    $N_occ      = undef;
	    $CA_occ     = undef;
	    $C_occ      = undef;
	    $O_occ      = undef;
	    $bb_defined = undef;
	}
    }
}
&abort ("no backbone in pdb file '$pdbfile'") if (! $bb_defined_somewhere);


# output
#
$outbuf = join ("\n", @out)."\n";
if ($outfile) {
#    print "creating $outfile\n";
    open (OUTFILE, '>'.$outfile);
    select (OUTFILE);
}
print $outbuf;
if ($outfile) {
    close (OUTFILE);
    select (STDOUT);
}
&abort ("failed to create out file '$outfile'") if ($outfile && ! -s $outfile);


# exit
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
    local $usage = qq{
usage: $0 
\t-pdbfile          <pdbfile>
\t[-outfile         <outfile>]
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "pdbfile=s", "outfile=s");


    # Check for legal invocation
    #
    if (! defined $opts{pdbfile}) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExistence ('f', $opts{pdbfile});	

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
    my $msg = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print STDERR "[$date]:$0:ABORT: $msg\n";
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
