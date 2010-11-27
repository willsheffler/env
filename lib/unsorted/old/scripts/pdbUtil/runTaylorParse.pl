#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.2 $
##  $Date: 2005/06/11 00:57:37 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$debug = 1;

$| = 1;                                              # disable stdout buffering

$resmap_dir     = $ENV{'BAKER_HOME'}."/dat/resmap";
$pdb_chain_dir  = $ENV{'BAKER_HOME'}."/dat/pdbs/pdb-chain";

$cleanPdb       = $ENV{'BAKER_HOME'}."/src/pdbUtil/cleanPdb.pl";

$taylor_bin_dir = $ENV{'BAKER_HOME'}."/src/shareware/taylor_domains/bin";

@taylor_bins    = ("dom", "best.csh", "record.csh");

@chains = qw (A B C D E F G H I J K L M N O P Q R S T U V W X Y Z 0 1 2 3 4 5 6 7 8 9);

###############################################################################
# init
###############################################################################

if ($#ARGV < 2) {
    print STDERR "usage: $0 <chain.ids> <bad_chains.ids> <out_dir>\n";
    exit -1;
}
$chains_file    = shift @ARGV;
$badchains_file = shift @ARGV;
$out_dir        = shift @ARGV;

if ($out_dir =~ /^\./) {
    $out_dir = `pwd`."/".$out_dir;
    $out_dir =~ s/\n//;
}
            
$chains_file_buf    = &bigFileBufArray ($chains_file);
$badchains_file_buf = &bigFileBufArray ($badchains_file);

###############################################################################
# main
###############################################################################

# READ BAD CHAINS
#
%badchain_ids = ();
for (my $id_i=0; $id_i <= $#{$badchains_file_buf}; ++$id_i) {
    ($chain_id, @rest) = split (/\s+/, $badchains_file_buf->[$id_i]);
    $badchain_ids{$chain_id} = 1;
}


# PROCESS CHAINS
#
for (my $id_i=0; $id_i <= $#{$chains_file_buf}; ++$id_i) {
    ($chain_id, @rest) = split (/\s+/, $chains_file_buf->[$id_i]);

    next if ($badchain_ids{$chain_id});

    print "\nDOING $chain_id\n";

    $chain_id =~ /^(\w)(\w\w)(\w)(\w)$/;
    $base = $1.$2.$3;
    $chain = $4;
    $folder = $2;

    # prep
    #
    $out_folder = "$out_dir/$folder/$chain_id";
    &runCmd (qq{rm -rf $out_folder})         if (-d $out_folder);
    &runCmd (qq{mkdir -p $out_folder})       if (! -d $out_folder);
    chdir $out_folder;


    # get clean pdb
    #
    &runCmd (qq{$cleanPdb -pdbfile $pdb_chain_dir/$folder/$chain_id.pdb.gz -outfile $out_folder/$chain_id.pdb}, "NODIE");


    # only do those with backbone density
    #
    if (! -f "$out_folder/$chain_id.pdb") {
        print "NO BACKBONE FOR chain_id\n";
        next;
    }
    else {

	# get resmap
	@residues = ();
	@resmap_seq2coord = ();
	foreach $line (&fileBufArray ("$resmap_dir/$folder/$chain_id.resmap")) {
	    ($res, $seq_n, $coord_n) = split (/\s+/, $line);
	    $residues[$seq_n]         = $res;
	    $resmap_seq2coord[$seq_n] = $coord_n;
	}

    
	# make bins local
	foreach $bin (@taylor_bins) {
	    system (qq{cp $taylor_bin_dir/$bin ./});
	}


	# run taylors (doesn't exit 0)
	&runCmd (qq{best.csh $chain_id.pdb 0 > taylor.out 2>&1}, 'NODIE');


	# tidy up UNK domain break records, and reorder domains sequentially
	if (-f "dom1.out") {
	    my @first_res = ();
	    for ($dom_n_reported=1; -f "dom$dom_n_reported.out"; ++$dom_n_reported) {
		foreach $line (&fileBufArray ("dom$dom_n_reported.out")) {
		    next if ($line =~ /UNK/);
		    next if ($line !~ /^ATOM/);
		    next if (substr ($line, 12, 4) ne ' CA ');
		    $res_n = substr ($line, 22, 4);
		    $res_n =~ s/\s+//g;
		    $first_res[$dom_n_reported-1] = $res_n;
		    last;
		}
	    }
	    my $dom_n_reported_sorted = &insertSortIndexList (\@first_res, 'increasing');
	    for (my $dom_i=0; $dom_i <= $#{$dom_n_reported_sorted}; ++$dom_i) {
		$dom_n = $dom_i + 1;
		$dom_n_reported = $dom_n_reported_sorted->[$dom_i] + 1;
		system (qq{fgrep -v UNK dom$dom_n_reported.out > dom$dom_n.out.clean});
	    }
	}


	# determine domains
	@domain_info = &getDomains ($chain_id, \@residues, \@resmap_seq2coord);


	# make taylor.pdb
	if (-f "dom1.out.clean") {
	    system (qq{cat dom*.out.clean > taylor.pdb});
	}
	else {
	    system (qq{grep ^ATOM $base$chain.pdb | fgrep ' CA ' > taylor.pdb});
	}


	# output determined domains
	open (DOMS, '>'."$chain_id.doms");
	print DOMS join ("\t", 'DOMAIN', 'COORDS', 'SEQ_RANGE') ."\n";
	for ($dom_n=1; $dom_n <= $#domain_info; ++$dom_n) {
	    $dom_id = ($#domain_info == 1) ? '_' : $dom_n;
	    print DOMS "$base$chain$dom_id\t$domain_info[$dom_n]\n";
	}
	close (DOMS);


	# clean up
	foreach $bin (@taylor_bins) {
	    system (qq{rm -f $bin});
	}
	system (qq{rm -rf doms});
	system (qq{rm -f topdom.log})  if (-f "topdom.log");
	system (qq{rm -f dom*.log})    if (-f "dom14.log");
	system (qq{rm -f dom*.out*})   if (-f "dom0.out");
    }
}

# done
exit 0;

###############################################################################
# subs
###############################################################################

sub getDomains {
    ($chain_id, $residues, $resmap_seq2coord) = @_;
    my @domain_info = ();
    my $line = '';
    my $i, $seg_i, $dom_n, $res_n;
    my $start_n, $stop_n;
    my $res_by_dom = +[];

    my $one_dom  = undef;
    my $no_parse = undef;
    my $log_id   = undef;
    foreach $line (&fileBufArray ('taylor.out')) {
	if ($line =~ /^ONE domain/) {
	    $one_dom = 'true';
	    last;
	}
	elsif ($line =~ /^NO match for $chain_id.pdb \( \d+ res.\) best at (\d+)/) {
	    $no_parse = 'true';
	    $log_id = $1;
	    $log_id = '10'  if ($log_id eq '0');
	    last;
	}
    }
    

    # find which residues are defined in each domain
    #
    if ($one_dom) {
        foreach $line (&fileBufArray ("$chain_id.pdb")) {
            next if ($line !~ /^ATOM/);
            next if (substr ($line, 12, 4) ne ' CA ');
            $res_n = substr ($line, 22, 4);
            $res_n =~ s/\s+//g;
	    # make mask
	    $res_by_dom->[1]->[$res_n] = 1;
	}
    }
    elsif ($no_parse) {
	foreach $line (&fileBufArray ("dom$log_id.log")) {
            if ($line =~ /^\d+\s*domains/) {
                $res_by_dom = +[];
                @domain_order = ();
            }
            if ($line =~ /^segment \d+ in domain (\d+) = (\d+) \(([A-Z])(\d+)\) --- (\d+) \(([A-Z])(\d+)\)/) {
                $dom_n_reported    = $1;
                $order_start_n     = $2;
                $coord_start_res   = $3;
                $coord_start_n     = $4;
                $order_stop_n      = $5;
                $coord_stop_res    = $6;
                $coord_stop_n      = $7;
		
                # goof balls don't number domains sequentially!
                $dom_n = undef;
                for ($dom_i=0; $dom_i <= $#domain_order; ++$dom_i) {
                    if ($domain_order[$dom_i] == $dom_n_reported) {
                        $dom_n = $dom_i + 1;
                        last;
                    }
                }
                if (! defined $dom_n) {
                    push (@domain_order, $dom_n_reported);
                    $dom_n = $#domain_order + 1;
                }

                # make mask
                for ($res_n=$coord_start_n; $res_n <= $coord_stop_n; ++$res_n) {
                    $res_by_dom->[$dom_n]->[$res_n] = 1;
                }
            }
	}
    
	# make dom?.out.clean files
	for ($dom_n=1; $dom_n <= $#{$res_by_dom}; ++$dom_n) {
	    my $out_buf = '';
	    foreach $line (&fileBufArray ("$chain_id.pdb")) {
		next if ($line !~ /^ATOM/);
		next if (substr ($line, 12, 4) ne ' CA ');
		$res_n = substr ($line, 22, 4);
		$res_n =~ s/\s+//g;
		if ($res_by_dom->[$dom_n]->[$res_n]) {
		    substr ($line, 21, 1) = $chains[$dom_n-1];
		    $out_buf .= "$line\n";
		}
	    }
	    open (OUT, '>'."dom$dom_n.out.clean");
	    print OUT $out_buf;
	    close (OUT);
	}
    }
    else {
	@domain_order = ();
	for ($dom_n=1; -f "dom$dom_n.out.clean"; ++$dom_n) {
	    foreach $line (&fileBufArray ("dom$dom_n.out.clean")) {
		next if ($line !~ /^ATOM/);
		next if (substr ($line, 12, 4) ne ' CA ');
		$res_n = substr ($line, 22, 4);
		$res_n =~ s/\s+//g;
		# make mask
		$res_by_dom->[$dom_n]->[$res_n] = 1;
	    }
	}
    }


    # add insertion mask
    #
    for ($dom_n=1; $dom_n <= $#{$res_by_dom}; ++$dom_n) {
        for ($i=1; $i <= $#{$res_by_dom->[$dom_n]}; ++$i) {
            if ($res_by_dom->[$dom_n]->[$i] == 1) {
                for ($dom_m=1; $dom_m <= $#{$res_by_dom}; ++$dom_m) {
                    next if ($dom_m == $dom_n);
                    $res_by_dom->[$dom_m]->[$i] = -1;
                }
            }
        }
    }
 
    # build domain_info
    #
    for ($dom_n=1; $dom_n <= $#{$res_by_dom}; ++$dom_n) {
        my @starts = ();
        my @stops  = ();
        $in_seg = undef;
        for ($res_n=1; $res_n <= $#{$res_by_dom->[$dom_n]}; ++$res_n) {
            next if ($res_by_dom->[$dom_n]->[$res_n] != 1);

            if (! $in_seg && $res_by_dom->[$dom_n]->[$res_n-1] != 1) {
                push (@starts, $res_n);
                $in_seg = 'true';
            }
            if ($in_seg && $res_by_dom->[$dom_n]->[$res_n+1] != 1) {
                $true_stop = 'true';
                for ($res_m=$res_n+1; $res_m <= $#{$res_by_dom->[$dom_n]}; ++$res_m) {
                    if ($res_by_dom->[$dom_n]->[$res_m] == 1) {
                        $true_stop = undef;
                        last;
                    }
                    elsif ($res_by_dom->[$dom_n]->[$res_m] == -1) {
                        last;
                    }
                }
                if ($true_stop) {
                    push (@stops, $res_n);
                    $in_seg = undef;
                }
            }
        }


        # build ranges now that we have them
        #
	my $coord_range = '';
        my $seq_range   = '';
        for ($seg_i=0; $seg_i <= $#starts; ++$seg_i) {
            $seq_start   = $starts[$seg_i];
            $seq_stop    = $stops[$seg_i];
            $res_start   = $residues->[$seq_start];
            $res_stop    = $residues->[$seq_stop];
            $coord_start = $resmap_seq2coord->[$seq_start];
            $coord_stop  = $resmap_seq2coord->[$seq_stop];

            $coord_range .= ($seg_i == 0) ? "$coord_start:$coord_stop"
                                          : ','."$coord_start:$coord_stop";
            $seq_range   .= ($seg_i == 0) ? "$res_start$seq_start:$res_stop$seq_stop"
                                          : ','."$res_start$seq_start:$res_stop$seq_stop";
        }

        $domain_info[$dom_n] = "$coord_range\t$seq_range";
    }                          


    return @domain_info;
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
    if ($file =~ /^\.Z$/ || $file =~ /\.gz$/) {
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
    if ($file !~ /^\.Z$/ && $file !~ /\.gz$/) {
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

# abort()
#
sub abort {
    my $msg = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print STDERR "[$date]:$0:ABORT: $msg\n";
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
    if ($file =~ /\.gz$|\.Z$/) {
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
    if ($file =~ /\.gz$|\.Z$/) {
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
    if ($file =~ /\.gz$|\.Z$/) {
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
