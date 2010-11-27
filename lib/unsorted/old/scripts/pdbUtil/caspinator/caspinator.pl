#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.6 $
##  $Date: 2004/07/14 01:57:26 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering

$man_group_id  = '4793-7654-3148';  # MANUAL: BAKER (CASP-6)
#$auto_group_id = '6750-1618-6576';  # AUTO:   BAKER-ROBETTA (CASP-6)
#$auto_group_id = 'BAKER-ROBETTA';  # AUTO:   BAKER-ROBETTA
$auto_group_id = '2596-3328-7939';  # AUTO:   BAKER-ROBETTA_04 (CASP-6)

$extra_header  = qq{
EXPDTA    THEORETICAL MODEL
AUTHOR    D.CHIVIAN,D.E.KIM,L.MALMSTROM,P.BRADLEY,T.ROBERTSON,P.MURPHY,
AUTHOR   2 C.E.M.STRAUSS,R.BONNEAU,C.A.ROHL,D.BAKER
JRNL        AUTH   D.CHIVIAN,D.E.KIM,L.MALMSTROM,P.BRADLEY,T.ROBERTSON,
JRNL        AUTH 2 P.MURPHY,C.E.M.STRAUSS,R.BONNEAU,C.A.ROHL,D.BAKER
JRNL        TITL   AUTOMATED PREDICTION OF CASP-5 STRUCTURES USING THE
JRNL        TITL 2 ROBETTA SERVER
JRNL        REF    PROTEINS: STRUCT., FUNCT.,    V.  53   524 2003
JRNL        REF  2 GENET.  
};
$extra_header =~ s/^\s+|\s+$//g;


$man_method_str  = 
qq{
ROSETTA provides both ab initio and
comparative models of protein domains. It
uses the ROSETTA fragment insertion method
[Simons et al. J Mol Biol 1997;268:209-225].
Comparative models are built from structures
detected by PSI-BLAST, FFAS03, or 3DJury-A1
and aligned by the K*SYNC alignment method.
Loop regions are assembled from fragments and 
optimized to fit the aligned template structure.
};

$auto_method_str =
qq{
ROBETTA provides both ab initio and
comparative models of protein domains. It
uses the ROSETTA fragment insertion method
[Simons et al. J Mol Biol 1997;268:209-225].
Comparative models are built from structures
detected by PSI-BLAST, FFAS03, or 3DJury-A1
and aligned by the K*SYNC alignment method.
Loop regions are assembled from fragments and 
optimized to fit the aligned template structure.
The procedure is fully automated.
------
ROBETTA_04 examines ensembles of alignments produced
parametrically with the K*SYNC alignment method with
multiple parents.  PSIBLAST level targets have
frozen templates plus loops modeled by fragments,
with models selected from the ensemble by full
heavy atom representation energetics.  FFAS03 and
3DJury level targets are allowed backbone flexibility
along the entire chain, including template regions,
with models selected from the ensemble by side-chain
centroid represention energetics.
};

#$res_shift = 0;
$x_shift   = 100;

###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
$conf             = $opts{conf};
$target_name      = $opts{targetname};
$model_num        = $opts{modelnum};
$cutdefs_file     = $opts{cutdefs};
$domainfiles_list = $opts{domainfileslist};
$bfactor_source   = $opts{bfactorsource};
$parents_list     = $opts{parentslist};
$keep_hetero      = $opts{keephetero};
$out_file         = $opts{outfile};

$outbuf = '';
@domainfiles = split (/,/, $domainfiles_list);

###############################################################################
# main
###############################################################################

# get group id
$group_id = ($conf =~ /^M/i) 
              ? $man_group_id
              : $auto_group_id;


# get method desc
$method_str = ($conf =~ /^M/i) 
              ? $man_method_str
              : $auto_method_str;
$method_str =~ s/^\s+|\s+$//g;


# clean target name
$target_name =~ s/^t(\d\d\d\d)$/T$1/;
$target_name =~ s/^t(\d\d\d)_?$/T0$1/i;


# read cutdefs
@cutdefs = &fileBufArray ($cutdefs_file);
for ($line_i=1; $line_i <= $#cutdefs; ++$line_i) {
    $domain_i = $line_i - 1;
    $cutdefs[$line_i] =~ s!^\s+|\s+$!!g;
    ($q_beg[$domain_i], 
     $q_end[$domain_i], 
     $q_len[$domain_i], 
     $m_beg[$domain_i], 
     $m_end[$domain_i], 
     $m_len[$domain_i], 
     $p_beg[$domain_i],
     $p_end[$domain_i],
     $p_id[$domain_i],
     $confidence[$domain_i],
     $p_src[$domain_i],
     $seq[$domain_i]) = split (/\s+/, $cutdefs[$line_i]);

     $parents[$domain_i] = (! $parents_list && 
			    $p_id[$domain_i] !~ /^n\/?a$/i &&
			    $p_src[$domain_i] !~ /^pfam$/i)
                             ? $p_id[$domain_i]
                             : 'N/A';

     # until rich's confidence function is here
     $confs[$domain_i] = (! $parents_list && 
			  $p_id[$domain_i] !~ /^n\/?a$/i &&
			  $p_src[$domain_i] !~ /^pfam$/i)
                           ? sprintf ("%-6.2f", $confidence[$domain_i])
                           : '0.00';

    # retain pfam family in p_src 
    if ($p_src[$domain_i] =~ /^pfam$/) {
	$p_src[$domain_i] .= ":$p_id[$domain_i]";
    }
}
# define linkers
for ($domain_i=1; $domain_i <= $#m_beg; ++$domain_i) {
    for ($res_n=$m_beg[$domain_i]; $res_n <= $m_end[$domain_i-1]; ++$res_n) {
	$linker[$res_n] = 'true';
    }
}


# read any bfactor source file
if ($bfactor_source) {
    foreach $line (&fileBufArray ($bfactor_source)) {
	next if ($line !~ /^ATOM/);
#	$atom_n = substr ($line, 6, 5);
#	$atom_n =~ s/\s+//g;
#	$bfactor_byatom[$atom_n] = substr ($line, 60, 6);
	$res_n  = substr ($line, 22, 5);
	$res_n  =~ s/\s+//g;
	$bfactor_byres[$res_n] = substr ($line, 60, 6);
    }
}


# get parents
if ($parents_list) {
    @cmdline_parents = split (/,/, $parents_list);
    for ($domain_i=0; $domain_i <= $#cmdline_parents; ++$domain_i) {
	if ($cmdline_parents[$domain_i]) {
	    @these_parents = split (/:/, $cmdline_parents[$domain_i]);
	    for ($parent_i=0; $parent_i <= $#these_parents; ++$parent_i) {
		$these_parents[$parent_i] =~ s/^(\w\w\w\w)(\w).*$/$1.'_'.$2/e;
		$these_parents[$parent_i] =~ s/__$//;
	    }
	    $parents[$domain_i] = join (' ', @these_parents);
	} else {
	    $parents[$domain_i] = 'N/A';
	}
    }
}


# start outbuf
$outbuf = qq{PFRMAT TS\nTARGET $target_name\nAUTHOR $group_id\n};
#$outbuf = qq{PFRMAT    TS\nTARGET    $target_name\n$extra_header\n};
foreach $line (split (/\n/, $method_str)) {
    $outbuf .= qq{METHOD    $line\n};
}            
$outbuf .= qq{MODEL     $model_num\n};


# zero hetatms
@hetatms = ();


# assemble pdbs
$atom_n = 1;

if ($#q_beg > 0 && $#domainfiles > 0) {           # haven't yet assembled pdbs
    
    for ($domain_i=0; $domain_i <= $#q_beg; ++$domain_i) {
	if ($domain_i == 0) {
	    @reported_parents = ();
	    for ($domain_j=0; $domain_j <= $#parents; ++$domain_j) {
		next if ($parents[$domain_j] eq 'N/A' && $#parents > 0);
		push (@reported_parents, $parents[$domain_j]);
	    }
	    push (@reported_parents, 'N/A')  if ($#reported_parents < 0);
	    $outbuf .= qq{PARENT }. join (", ", @reported_parents) .qq{\n};  
	}
	$outbuf .= qq{REMARK PARENT    $parents[$domain_i]\n};  
	$outbuf .= qq{REMARK PARSRC    $p_src[$domain_i]\n};  
	$outbuf .= qq{REMARK SCORE     $confs[$domain_i]\n};  
	
	$global_shift = $q_beg[$domain_i];
	$local_shift  = $q_beg[$domain_i] - $m_beg[$domain_i] + 1; 
	$start_res    = $q_beg[$domain_i] - $m_beg[$domain_i] + 1;
	$stop_res     = $q_end[$domain_i] - $m_beg[$domain_i] + 1; 
	
	# body
	$pdb_buf = &bigFileBufArray ($domainfiles[$domain_i]);
	for ($i=0; $i <= $#{$pdb_buf}; ++$i) {
	    if ($pdb_buf->[$i] =~ /^ATOM/) {
		
		# fix any missing occupancy and temp
		$pdb_buf->[$i] .= ' 'x(80 - length($pdb_buf->[$i]));
		$occ  = substr ($pdb_buf->[$i], 54, 6);
		$temp = substr ($pdb_buf->[$i], 60, 6);
		substr ($pdb_buf->[$i], 54, 6) = sprintf ("%6.2f", 1.00)  if ($occ =~ /^\s*$/ || $occ <= 0);
		substr ($pdb_buf->[$i], 60, 6) = sprintf ("%6.2f", 10.00) if ($temp =~ /^\s*$/ || $parents[$domain_i] eq 'N/A');
		
		# adjust residue numbering if necessary, and skip extra coordinates
		$resnum = substr ($pdb_buf->[$i], 22, 5);
		next if ($resnum < $start_res || $resnum > $stop_res);
		
		# fix residue and atom numbering
		$res_n = $resnum - $local_shift + $global_shift;
		substr ($pdb_buf->[$i], 22, 4) = sprintf ("%4d", $res_n);
		substr ($pdb_buf->[$i], 5, 6)  = sprintf ("%6d", $atom_n++);
		
		# color linkers
		if ($linker[$res_n]) {
		    $temp = substr ($pdb_buf->[$i], 60, 6);
		    substr ($pdb_buf->[$i], 60, 6) = sprintf ("%6.2f", $temp+5.0);
		}

		# adjust x coordinate
		$x = substr ($pdb_buf->[$i], 30, 8);
		$x += $x_shift * $domain_i;
		substr ($pdb_buf->[$i], 30, 8) = sprintf ("%8.3f", $x);
		
		# append
		$outbuf .= $pdb_buf->[$i]."\n";
	    }
	    elsif ($pdb_buf->[$i] =~ /^HETATM/) {
		substr ($pdb_buf->[$i], 60, 6) = sprintf ("%6.2f", 5.00); 
		push (@hetatms, $pdb_buf->[$i]);
	    }
	}
#	$outbuf .= qq{TER\n};
    }
}
else {                                           # using a single assembled pdb
	
    # body
    $pdb_buf = &bigFileBufArray ($domainfiles[0]);
    for ($i=0; $i <= $#{$pdb_buf}; ++$i) {
	if ($pdb_buf->[$i] =~ /^ATOM/) {
	    
	    # fix any missing occupancy and temp
	    $pdb_buf->[$i] .= ' 'x(80 - length($pdb_buf->[$i]));
	    $occ  = substr ($pdb_buf->[$i], 54, 6);
	    $temp = substr ($pdb_buf->[$i], 60, 6);
	    substr ($pdb_buf->[$i], 54, 6) = sprintf ("%6.2f", 1.00)   if ($occ =~ /^\s*$/ || $occ <= 0);
	    substr ($pdb_buf->[$i], 60, 6) = sprintf ("%6.2f", 10.00)  if ($temp =~ /^\s*$/);
	    
	    # find which domain we're in, and print domain_header if new
	    $resnum = substr ($pdb_buf->[$i], 22, 5);
	    for ($domain_i=0; $domain_i <= $#q_beg; ++$domain_i) {
		if ($resnum == $q_beg[$domain_i] && ! $header_printed[$domain_i]) {
		    $header_printed[$domain_i] = 'true';
#		    if ($resnum != 1) {
#			$outbuf .= qq{TER\n};
#		    }
		    if ($domain_i == 0) {
			@reported_parents = ();
			for ($domain_j=0; $domain_j <= $#parents; ++$domain_j) {
			    next if ($parents[$domain_j] eq 'N/A' && $#parents > 0);
			    push (@reported_parents, $parents[$domain_j]);
			}
			push (@reported_parents, 'N/A')  if ($#reported_parents < 0);
			$outbuf .= qq{PARENT }. join (", ", @reported_parents) .qq{\n};  
		    }
		    $outbuf .= qq{REMARK PARENT    $parents[$domain_i]\n};
		    $outbuf .= qq{REMARK PARSRC    $p_src[$domain_i]\n};  
		    $outbuf .= qq{REMARK SCORE     $confs[$domain_i]\n};  
		}
	    }
		
	    # restore bfactor if there's a source
#	    substr ($pdb_buf->[$i], 60, 6) = sprintf ("%6.2f", $bfactor_byatom[$atom_n])  if ($bfactor_source);
	    substr ($pdb_buf->[$i], 60, 6) = sprintf ("%6.2f", $bfactor_byres[$resnum])  if ($bfactor_source);
		
	    # fix atom numbering
	    substr ($pdb_buf->[$i], 5, 6)  = sprintf ("%6d", $atom_n++);
		
            # append
	    $outbuf .= $pdb_buf->[$i]."\n";
	}
	elsif ($pdb_buf->[$i] =~ /^HETATM/) {
	    substr ($pdb_buf->[$i], 60, 6) = sprintf ("%6.2f", 5.00); 
	    push (@hetatms, $pdb_buf->[$i]);
	}
    }
#    $outbuf .= qq{TER\n};
}


# chain is done
#
$outbuf .= sprintf ("TER  %6d\n", $atom_n++);


# renumber HETATMS and add them too
#
if ($keep_hetero) {
    for ($i=0; $i <= $#hetatms; ++$i) {
	substr ($hetatms->[$i], 5, 6) = sprintf ("%6d", $atom_n++);
	$outbuf .= $hetatms[$i]."\n";
    }
}


# file is done
#
$outbuf .= qq{END\n};


# ADD:
# rasmol HELIX, SHEET, TURN (e.g. charlie's add_SS.pl)?
# energies?
# torsion angles?


# output
if ($out_file) {
    open (OUT, '>'.$out_file);
    select (OUT);
}
print $outbuf;
if ($out_file) {
    close (OUT);
    select (STDOUT);
}

# done
exit 0;

###############################################################################
# subs
###############################################################################

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{usage: $0
\t -targetname      <target_name>
\t -modelnum        <model_num>
\t -cutdefs         <cutdefs_file>
\t -domainfileslist <domainfile1,domainfile2,...>
\t[-conf            <MANUAL/AUTOMATIC>]               (def: MANUAL)
\t[-bfactorsource   <bfactor_source_pdb>]             (def: none)
\t[-parentslist     <parent1,parent2a:parent2b,...>]  (def: reads cutdefs)
\t[-keephetero      <TRUE/FALSE>]                     (def: FALSE)
\t[-outfile         <out_file>]                       (def: STDOUT)
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "targetname=s",
		 "modelnum=i",
		 "cutdefs=s",
		 "domainfileslist=s",
		 "conf=s",
		 "bfactorsource=s",
		 "parentslist=s",
		 "keephetero=s",
		 "outfile=s"
		 );

    # Check for legal invocation
    if (! defined $opts{targetname} ||
	! defined $opts{modelnum} ||
	! defined $opts{cutdefs} ||
	! defined $opts{domainfileslist}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }

    # Defaults
    $opts{conf}       = 'MANUAL'  if (! defined $opts{conf});
    $opts{keephetero} = undef     if ($opts{keephetero} =~ /^F/i);

    # Existence checks
    &checkExist ('f', $opts{cutdefs});
    &checkExist ('f', $opts{bfactorsource})  if (defined $opts{bfactorsource});

    return %opts;
}

###############################################################################
# util
###############################################################################

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
    my ($dir, $fullpath_flag) = shift;
    my $inode;
    my @inodes = ();
    my @files = ();
    
    opendir (DIR, $dir);
    @inodes = readdir (DIR);
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
    my $cmd = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print "[$date]:$0: $cmd\n" if ($debug);
    &abort ("failure running cmd: $cmd\n") if (system ($cmd) != 0);
    #if (system ($cmd) != 0) {
    #    print STDERR ("failure running cmd: $cmd\n");
    #    return -2;
    #}
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
    print STDERR "[$date]:$0: $msg\n";
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
