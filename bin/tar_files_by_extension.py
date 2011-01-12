#!/usr/bin/env python
import sys,os,subprocess

EXCLUDED_EXT = ("tar","tgz","gz","bz2","zip","rar")

def files_with_ext(ext,names):
	r = []
	if ext == "NOEXTENSION":
		r = [n for n in names if n.count(".") == 0]
	else:
		x = "."+ext
		for n in names:
			if not n.endswith(x): continue
			if n == x: continue
			r.append(n)
	return r
	
def split_files(files,N):
	batches = []
	for i in range(999999):
		x = files[N*i:N*(i+1)]
		if x: batches.append( x )
		else: break
	return batches

def execute_and_check(cmd,cwd=None):
	if len(cmd) > 200: print "tar_files_by_extension.py:",cmd[:99],"...",cmd[-99:]
	else:              print "tar_files_by_extension.py:",cmd
	p = subprocess.Popen(cmd, cwd=cwd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	so = p.stdout.read()
	se = p.stderr.read()
	if so or se:
		print "executing cmd '"+cmd+"' failed!" 
		print "stdout:",so
		print "stderr:",se
		sys.exit()

def sanitize_fname(fn):
	return fn#.replace("'",r"\'").replace("(",r"\(").replace(")",r"\)")

def filelist_nchar(files):
	return len("'"+"' '".join(files)+"'")

def get_extensions(names):
	ext = set()
	for n in names:
		if n.count(".") == 0: continue # ignore if no extension
		if n.startswith("."): continue # no dotfiles
		if n.endswith(".gz"):
			if n.count(".") == 1: continue # don't bother with abcd.gz
			x = ".".join( n.split(".")[-2:] )
			if n != "."+x: ext.add( x )
			# print n,x
		else:
			x = n.split(".")[-1]
			if n != "."+x: ext.add( x )
			# print n,x
	ext = set([x for x in ext if not x in EXCLUDED_EXT])
	ext.add("NOEXTENSION")
	return ext

def visit(arg,dirname,names_in):
	print
	print dirname
	names = [n for n in names_in if not os.path.isdir(dirname+"/"+n)]
	exts = get_extensions(names)
	print exts
	for ext in exts:
		if ext in EXCLUDED_EXT: continue
		files = files_with_ext(ext,names)
		if not files: continue	
		N = 30000
		batches = split_files(files,N)
		while True:
			m = max([filelist_nchar(b) for b in batches])
			if m > 250000: 
				print "reducing batch size to",N
				N -= 1000
				batches = split_files(files,N)				
			else:
				break
		newext = ext.replace(".","_")	
		print "BATCH",ext,len(batches),[len(x) for x in batches]
		for i,batch in enumerate(batches):
			pref = "_%03i."%i
			if len(batches) == 1: pref = ""
			fnames = '"'+'" "'.join([sanitize_fname(x) for x in batch])+'"'
			execute_and_check( "tar -f _TAR_BY_EXTENSION%s.%s.tgz -zc %s"%(pref,newext,fnames) , dirname )
			execute_and_check( "rm  -f                                %s"%(            fnames) , dirname )
		# else:
		# 	filebatch = split_files(files)
		# 	execute_and_check(                     "tar -f %s/%s.tar -c %s/%s"%(dirname,ext,dirname,files[0]) )
		# 	for n in files[1:]:	execute_and_check( "tar -f %s/%s.tar -r %s/%s"%(dirname,ext,dirname,    n   ) )
		# 	execute_and_check(                     "gzip   %s/%s.tar"         %(dirname,ext                 ) )
		# 	for n in files:     execute_and_check( "rm -f               %s/%s"%(            dirname,    n   ) )
		# 	execute_and_check("mv %s/%s.tar.gz %s/%s.tgz"%(dirname,ext,dirname,ext))

def main():
	for path in sys.argv[1:]:
		os.path.walk(path,visit,None)

if __name__ == '__main__':
	main()