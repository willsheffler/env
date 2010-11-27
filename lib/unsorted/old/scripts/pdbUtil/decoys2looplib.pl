#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.5 $
##  $Date: 2004/07/27 06:06:13 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering

$max_splice_thresh = 101;

###############################################################################
# init
###############################################################################

# argv
#my %opts = &getCommandLineOptions ();
#my $pdbfile = $opts{pdbfile};
#my $code4   = $opts{code4};
#my $outfile = $opts{outfile};

if ($#ARGV < 2) {
#    print STDERR "usage: $0 <zone_file> <splice_cutoff> <decoys...>\n";
    print STDERR "usage: $0 <old_looplib> <splice_cutoff> <decoys...>\n";
    exit -1;
}
#$zone_file        = shift @ARGV;
$old_looplib_file = shift @ARGV;
$splice_cutoff    = shift @ARGV;
@decoys = @ARGV;

$splice_cutoff_incr = $splice_cutoff / 2.0;

###############################################################################
# main
###############################################################################

#
# NEVER USE DECOYS THAT WERE CREATED WITH DIFFERENT ZONE FILE!!!
#

# determine query length
#
foreach $line (&fileBufArray ($decoys[0])) {
    if ($line =~ /^ATOM/) {
	$res_n = substr ($line, 22, 4);
    }
}
$seq_len = $res_n;


# read loop definitions from zones (don't need with old_looplib)
#
#@zones    = &fileBufArray ($zone_file);
#@loopdefs = &zones2loops ($seq_len, @zones);
#foreach $loopdef (@loopdefs) {
#    $loopdef =~ s/\s+//g;
#    if ($loopdef =~ /^loop(\d+)-(\d+):(\w+)-(\w+)$/i) {
#	$q_loop_start = $1;
#	$q_loop_stop = $2;
#	push (@loop_starts, $q_loop_start);
#	push (@loop_stops,  $q_loop_stop);
#    }	
#}


# read loop definitions from old_looplib
#
@loop_starts = ();
@loop_stops  = ();
@loop_dirs   = ();
$last_loop_start = undef;
foreach $loop_frag_line (&fileBufArray ($old_looplib_file)) {
    $loop_start = substr ($loop_frag_line, 4, 4);
    $loop_stop  = substr ($loop_frag_line, 9, 4);
    $loop_dir   = substr ($loop_frag_line, -3, 3);

    $loop_start =~ s/\s+//g;
    $loop_stop  =~ s/\s+//g;
    $loop_dir   =~ s/\s+//g;

    if (! defined $last_loop_start || $loop_start != $last_loop_start) {
	$last_loop_start = $loop_start;
	push (@loop_starts, $loop_start);
	push (@loop_stops,  $loop_stop);
	push (@loop_dirs,   $loop_dir);
#	print STDERR "$loop_start\t$loop_stop\t$loop_dir\n";  # DEBUG
    }
}


# read decoys
#
for ($loop_i=0; $loop_i <= $#loop_starts; ++$loop_i) {
    $splice_cutoff_for_loop[$loop_i] = $splice_cutoff;
}
$all_loops_close = undef;
while ($splice_cutoff < $max_splice_thresh && ! $all_loops_close) {
    $looplib = +{};
    @closing_loop_found = ();
    $decoy_cnt = 0;
    foreach $decoy (@decoys) {
	
	++$decoy_cnt;
	print STDERR "$decoy_cnt decoys processed\n"  if (($decoy_cnt % 100) == 0);
	
	$pdb_file_buf = &bigFileBufArray ($decoy);
	
	@splicerms = ();
	@angle_buf = ();
	$started   = undef;
	
	# get splicerms
	foreach $line (@$pdb_file_buf) {
	    next if ($line !~ /splicerms/);
	    $line =~ s/^\s*splicerms\:\s*|\s+$//g;
	    @splicerms = split (/\s+/, $line);
	    last;
	}
	
	# get angles
	foreach $line (@$pdb_file_buf) {
	    if ($line =~ /^complete/) {
		$started = 'true';
		next;
	    }
	    next if (! $started);
	    push (@angle_buf, $line);
	    
	    last if ($line !~ /^\s*\d+ [HEL]\s+[\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+/);
	}
	
	# get loops
	$loop_info = +[];
	for ($loop_i=0; $loop_i <= $#loop_starts; ++$loop_i) {
	    $ss_list   = +[];
	    $phi_list  = +[];
	    $psi_list  = +[];
	    $omg_list  = +[];
	    $loop_start = $loop_starts[$loop_i];
	    $loop_stop  = $loop_stops[$loop_i];
	    $loop_dir   = $loop_dirs[$loop_i];
	    
	    for ($res_i=0; $loop_start+$res_i <= $loop_stop; ++$res_i) {
		$angle_buf[$loop_start-1+$res_i] =~ s/^\s+|\s+$//g;
		($res_n, $ss, $phi, $psi, $omg, $loop_src, $dev) = split (/\s+/, $angle_buf[$loop_start-1+$res_i]);
		
		push (@$ss_list,  $ss);
		push (@$phi_list, $phi);
		push (@$psi_list, $psi);
		push (@$omg_list, $omg);
	    }
	    
	    $loop_info->[$loop_i]->{splicerms}  = $splicerms[$loop_i];
	    $loop_info->[$loop_i]->{loop_start} = $loop_start;
	    $loop_info->[$loop_i]->{loop_stop}  = $loop_stop;
	    $loop_info->[$loop_i]->{loop_dir}   = $loop_dir;
	    $loop_info->[$loop_i]->{ss_list}    = $ss_list;
	    $loop_info->[$loop_i]->{phi_list}   = $phi_list;
	    $loop_info->[$loop_i]->{psi_list}   = $psi_list;
	    $loop_info->[$loop_i]->{omg_list}   = $omg_list;
	}
	
	# attach loops to looplib
	for ($loop_i=0; $loop_i <= $#{$loop_info}; ++$loop_i) {
	    $line = '';
	    $loop_start    = $loop_info->[$loop_i]->{loop_start};
	    $loop_stop     = $loop_info->[$loop_i]->{loop_stop};
	    $loop_dir      = $loop_info->[$loop_i]->{loop_dir};
	    $loop_len      = $loop_stop - $loop_start + 1;
	    $ss_str        = join ('', @{$loop_info->[$loop_i]->{ss_list}});
	    $splicerms_val = $loop_info->[$loop_i]->{splicerms};
	    
	    next if ($splicerms_val > $splice_cutoff_for_loop[$loop_i]);
	    $closing_loop_found[$loop_i] = 'true';
	    
	    $line = sprintf ("%3d %4d %4d %s %9.4f %9.4f %76s %4s ", $loop_len, $loop_start, $loop_stop, $ss_str, $splicerms_val, 0.0, '', ('0'x(4-length($decoy_cnt))).$decoy_cnt);
	    
	    for ($res_i=0; $res_i <= $#{$loop_info->[$loop_i]->{phi_list}}; ++$res_i) {
		$line .= sprintf ("%7.2f %7.2f %7.2f   ", $loop_info->[$loop_i]->{phi_list}->[$res_i], $loop_info->[$loop_i]->{psi_list}->[$res_i], $loop_info->[$loop_i]->{omg_list}->[$res_i]);
	    }
	    
	    if ($loop_dir eq '  1' || $loop_dir eq ' -1') {
		$line .= sprintf (" %2d", $loop_dir);
	    }
	    
	    push (@{$looplib->{$loop_start}}, $line);
	}
    }

    # ensure we have closers at each position
    $all_loops_close = 'true';
    for ($loop_i=0; $loop_i <= $#loop_starts; ++$loop_i) {
	if (! defined $closing_loop_found[$loop_i]) {
	    &alert ("no closing loops found at $splice_cutoff A for loop ".($loop_i+1)."... increasing threshold for this loop to ".($splice_cutoff+$splice_cutoff_incr));
	    $splice_cutoff_for_loop[$loop_i] += $splice_cutoff_incr;
	    $all_loops_close = undef;
	}
    }
    if (! $all_loops_close) {
	$splice_cutoff += $splice_cutoff_incr;
	next;
    } 
}


# write output
#
foreach $loop_start (sort {$a<=>$b} keys %{$looplib}) {
    print join ("\n", @{$looplib->{$loop_start}})."\n";
}


# done
exit 0;

###############################################################################
# subs
###############################################################################

# zones2loops
#
sub zones2loops {
    my ($seq_len, @zones) = @_;
    my @loops = ();
    my @q_starts = ();
    my @q_stops = ();
    my @p_starts = ();
    my @p_stops = ();
    my $line = undef;

    foreach $line (@zones) {
        $line =~ s/zone//;
        $line =~ s/\s+//g;
        ($q_start, $q_stop, $p_start, $p_stop) = split (/[^\d]+/, $line);
        push (@q_starts, $q_start);
        push (@q_stops, $q_stop);
        push (@p_starts, $p_start);
        push (@p_stops, $p_stop);
    }

    for (my $i=0; $i <= $#q_starts; ++$i) {
        if ($i == 0 && $q_starts[0] != 1) {
            $q_loop_start = 1;
            $q_loop_stop  = $q_starts[0]-1;
            $p_loop_start = 'NA';
            $p_loop_stop  = $p_starts[0]-1;                   
            push (@loops, "loop $q_loop_start-$q_loop_stop:$p_loop_start-$p_loop_stop");
        }
        if ($i > 0) {
            $q_loop_start = $q_stops[$i-1]+1;
            $q_loop_stop  = $q_starts[$i]-1;
            $p_loop_start = $p_stops[$i-1]+1;
            $p_loop_stop  = $p_starts[$i]-1;
            push (@loops, "loop $q_loop_start-$q_loop_stop:$p_loop_start-$p_loop_stop");
        }
        if ($i == $#q_starts && $q_stops[$i] != $seq_len) {
            $q_loop_start = $q_stops[$i]+1;
            $q_loop_stop  = $seq_len;
            $p_loop_start = $p_stops[$i]+1;
            $p_loop_stop  = 'NA';
            push (@loops, "loop $q_loop_start-$q_loop_stop:$p_loop_start-$p_loop_stop");
        }
    }

    return @loops;
} 


# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{usage: $0 -pdbfile <pdb_file> -code4 <code4> -outfile <looplib_file>\n};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, "pdbfile=s", "code4=s", "outfile=s");

    # Check for legal invocation
    if (! defined $opts{pdbfile} || 
	! defined $opts{code4}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{pdbfile});

    return %opts;
}

###############################################################################
# util
###############################################################################

# iterElimSortIndexList ()
#
sub iterElimSortIndexList {
    my ($val1_list, $val2_list, $fraction, $direction) = @_;
    my $index_list        = +[];
    my $local_index_list  = +[];
    my $local_val_list    = +[];
    my $local_sorted_list = +[];
    my ($index, $i, $j);

    my $sorted_val1_list = &insertSortIndexList ($val1_list, $direction);
    for ($i=0; $i <= $#{$sorted_val1_list}; ++$i) {
	$index_list->[$i] = $sorted_val1_list->[$i];
    }

    my $done = undef;
    my $toggle = 2;
    $cut = int ($#{$index_list} * $fraction);
    $last_cut = $#{$index_list};
    while ($cut > 0) {
	# sort the right half ("discards")
	$local_index_list = +[];
	$local_val_list   = +[];
	for ($j=0; $cut+$j+1 <= $last_cut; ++$j) {
	    $index                  = $index_list->[$cut+$j+1];
	    $local_index_list->[$j] = $index;
	    $local_val_list->[$j]   = ($toggle == 1) ? $val1_list->[$index]
		                                     : $val2_list->[$index];
	}
	$local_sorted_index_list = &insertSortIndexList ($local_val_list, $direction);
	for ($j=0; $cut+$j+1 <= $last_cut; ++$j) {
	    $local_index = $local_sorted_index_list->[$j];
	    $index_list->[$cut+$j+1] = $local_index_list->[$local_index];
	}

	# sort the left half ("keeps")
	$local_index_list = +[];
	$local_val_list   = +[];
	for ($j=0; $j <= $cut; ++$j) {
	    $index                  = $index_list->[$j];
	    $local_index_list->[$j] = $index;
	    $local_val_list->[$j]   = ($toggle == 1) ? $val1_list->[$index]
		                                     : $val2_list->[$index];
	}
	$local_sorted_index_list = &insertSortIndexList ($local_val_list, $direction);
	for ($j=0; $j <= $cut; ++$j) {
	    $local_index = $local_sorted_index_list->[$j];
	    $index_list->[$j] = $local_index_list->[$local_index];
	}
	
	# update cut and toggle
	$toggle = ($toggle == 1) ? 2 : 1;
	$last_cut = $cut;
	$cut = int ($last_cut * $fraction);
    }

    return $index_list;
}

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

# readFiles
#
sub readFiles {
    my ($dir, $fullpath_flag) = @_;
    my $inode;
    my @inodes = ();
    my @files = ();
    
    opendir (DIR, $dir);
    @inodes = sort readdir (DIR);
    closedir (DIR);
    foreach $inode (@inodes) {
	next if (! -f "$dir/$inode");
	next if ($inode =~ /^\./);
	push (@files, ($fullpath_flag) ? "$dir/$inode" : "$inode");
    }
    return @files;
}

# createDir
#
sub createDir {
    my $dir = shift;
    if (! -d $dir && (system (qq{mkdir -p $dir}) != 0)) {
	print STDERR "$0: unable to mkdir -p $dir\n";
	exit -2;
    }
    return $dir;
}

# copyFile
#
sub copyFile {
    my ($src, $dst) = @_;
    if (system (qq{cp $src $dst}) != 0) {
	print STDERR "$0: unable to cp $src $dst\n";
	exit -2;
    }
    return $dst;
}

# zip
#
sub zip {
    my $file = shift;
    if ($file =~ /^\.Z/ || $file =~ /\.gz/) {
	print STDERR "$0: ABORT: already a zipped file $file\n";
	exit -2;
    }
    if (system (qq{gzip -9 $file}) != 0) {
	print STDERR "$0: unable to gzip -9 $file\n";
	exit -2;
    }
    $file .= ".gz";
    return $file;
}

# unzip
#
sub unzip {
    my $file = shift;
    if ($file !~ /^\.Z/ && $file !~ /\.gz/) {
	print STDERR "$0: ABORT: not a zipped file $file\n";
	exit -2;
    }
    if (system (qq{gzip -d $file}) != 0) {
	print STDERR "$0: unable to gzip -d $file\n";
	exit -2;
    }
    $file =~ s/\.Z$|\.gz$//;
    return $file;
}

# remove
#
sub remove {
    my $inode = shift;
    if (system (qq{rm -rf $inode}) != 0) {
	print STDERR "$0: unable to rm -rf $inode\n";
	exit -2;
    }
    return $inode;
}
     
# runCmd
#
sub runCmd {
    my ($cmd, $nodie) = @_;
    my $ret;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print "[$date]:$0:RUNNING: $cmd\n" if ($debug);
    $ret = system ($cmd);
    #$ret = ($?>>8)-256;
    if ($ret != 0) {
	$date = `date +'%Y-%m-%d_%T'`;  chomp $date;
	print STDERR ("[$date]:$0: FAILURE (exit: $ret): $cmd\n");
	if ($nodie) {
	    return $ret;
	} else {
	    exit $ret;
	}
    }
    return 0;
}

# logMsg()
#
sub logMsg {
    my ($msg, $logfile) = @_;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    else {
	select (STDERR);
    }
    print "[$date]:$0: $msg\n";
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
