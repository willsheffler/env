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

###############################################################################
# init
###############################################################################

$| = 1;                                                   # don't buffer stdout

if ($#ARGV < 2) {
    print STDERR "usage: $0 <res_i> <res_j> <pdbfile1> [<pdbfile2>...]\n";
    exit -1;
}
$res_i    = shift @ARGV;
$res_j    = shift @ARGV;
@pdbfiles = @ARGV;

###############################################################################
# main
###############################################################################

# read
#
$cnt = 0;
foreach $pdbfile (@pdbfiles) {
    ++$cnt;
    print STDERR "reading $cnt"."th pdb\n"  if (($cnt % 10) == 0);
    @buf = &fileBufArray ($pdbfile);

    $res_i_found = undef;
    $res_j_found = undef;

    for ($i=0; $i <= $#buf; ++$i) {
	if ($buf[$i] =~ /^ATOM/) {
	    $atomtype = substr ($buf[$i], 13, 2);
	    $res_n    = substr ($buf[$i], 22, 4);
	    if ($atomtype eq 'CA') {
		if ($res_i == $res_n) {
		    $res_i_found = 'true';
		    $x_i = substr ($buf[$i], 30, 8);
		    $y_i = substr ($buf[$i], 38, 8);
		    $z_i = substr ($buf[$i], 46, 8);
		}
		elsif ($res_j == $res_n) {
		    $res_j_found = 'true';
		    $x_j = substr ($buf[$i], 30, 8);
		    $y_j = substr ($buf[$i], 38, 8);
		    $z_j = substr ($buf[$i], 46, 8);
		}
		if ($res_i_found && $res_j_found) {
		    $dist = sqrt (($x_i-$x_j)*($x_i-$x_j)+($y_i-$y_j)*($y_i-$y_j)+($z_i-$z_j)*($z_i-$z_j));
		    printf ("%s\t%8.3f\n", $pdbfile, $dist);
		    last;
		}
	    }
	}
    }
}

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
    local $usage = qq{usage: $0 -pdbfile <pdbfile> -chain <chain>};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "pdbfile=s", "chain=s");


    # Check for legal invocation
    #
    if (! defined $opts{pdbfile}) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExistence ('f', $opts{pdbfile});	


    # Defaults
    #
    $opts{chain} = ' '  if (! defined $opts{chain} || $opts{chain} eq '0' || $opts{chain} eq '_');
    $opts{chain} = uc $opts{chain};
 
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
