#!/usr/bin/perl
#----------- insert this section to exclude constraints on glycine -------
die ("Usage: $0 <constraint file>") if (@ARGV<1);  
$base_cst = "$ARGV[0]";
#$generation = "$ARGV[1]";
#$base_code = "$ARGV[2]";

$base_cst =~ /^(....)(.)/;
$id = $1;
$chain = $2;
$pdbfile_name = $id.".pdb";
#$base_cst = "1erv_-1t7pB.3dpair.set/1erv_-1t7pB.3dpair.base.cst";
open (IN_PDB, $pdbfile_name) || die "can't open pdb file!";
open (IN_CST, $base_cst) || die "can't open constraint file!";
open (OUT, ">tmp_cst.out") || die "can't open tmp cst file!";

%residue_pos = ();
@pdb = <IN_PDB>;
foreach $pdb (@pdb){
    if (($pdb =~ /^ATOM/) && ($pdb =~ /CA/)) {
	$resnum_tmp = substr ($pdb, 22, 4);   #get 4 spaces so some will be spaces
	$resnum_tmp =~ /\s+(\d+)/;            #extract the non-white spaces to make
	$resnum = $1;                         #sure only numbers are stored in hash
	$restype = substr ($pdb, 17, 3);
	$residue_pos{$resnum} = $restype;
    }
}


@cst = <IN_CST>;
foreach $cst (@cst) {
    if ($cst =~ /^\S+/) {
	print OUT $cst;
    }
    elsif ($cst =~ /\s+(\d+)\s+CB\s+(\d+)\s+CB/) {
	$res1 = $1;
	$res2 = $2;
	$type = $residue_pos{$res1};
	if (($residue_pos{$res1} !~ /GLY/) && ($residue_pos{$res2} !~ /GLY/)) {
	    print OUT $cst;
	}
    }
}

system "mv tmp_cst.out $base_cst";
