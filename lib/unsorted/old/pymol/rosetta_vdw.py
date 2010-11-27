from pymol import cmd


def useRosettaRadii():
	cmd.alter("element C", "vdw=2.00")
	cmd.alter("element N", "vdw=1.75")
	cmd.alter("element O", "vdw=1.55")
	cmd.alter("element H", "vdw=1.00")
	cmd.alter("element P", "vdw=1.90")
	cmd.set("sphere_scale", 1.0)
	

cmd.extend('useRosettaRadii', useRosettaRadii)
