if [[ $ZSH_VERSION == 4.2.6 ]]; then
	if [[ -e /usr/local/bin/zsh ]]; then
		echo "RUN NEW ZSH"
		/usr/local/bin/zsh;
		exit
	fi
fi

ZSH=$HOME/.oh-my-zsh
	kver=`uname -r`
	if [[ $kver[-4,-2] == '.el' ]]; then
		# no git
		# ZSH_THEME="cypher"
		# ZSH_THEME="dieter"
		ZSH_THEME="evan"
		# ZSH_THEME="imajes"
		plugins=(colorize colored-man zsh-syntax-highlighting cp dircycle pip)
		#function st { ssh $MYIP "cd `pwd|sed -e 's=/work/sheffler=/baker=g'`; '/Applications/Sublime Text.app/Contents/SharedSupport/bin/subl' $@" }
		MYIP=172.16.131.130
	else
		ZSH_THEME="wedisagree"
		#ZSH_THEME="agnoster"
		plugins=(brew colorize colored-man zsh-syntax-highlighting cp dircycle osx pip sublime) # git git-fast git-extras
	fi

source $ZSH/oh-my-zsh.sh
	export LC_ALL=en_US.UTF-8
	export LANG=en_US.UTF-8
	unsetopt correct_all
	_per_directory_history_is_global=true


export     PATH="."
	export PATH=$PATH:$HOME/local/bin
	export PATH=$PATH:$HOME/sw/bin
	export PATH=$PATH:$HOME/bin
	export PATH=$PATH:$HOME/lib/scripts
	export PATH=$PATH:$HOME/bin
	export PATH=$PATH:/usr/local/bin
	export PATH=$PATH:/usr/bin
	export PATH=$PATH:/bin
	export PATH=$PATH:/usr/sbin
	export PATH=$PATH:/sbin
	export PATH=$PATH:/usr/X11/bin
	export PATH=$PATH:$HOME/pymol
	export PATH=$PATH:$HOME/rosetta/matdes/source/cmake/build_omp


	if [[ $kver[-4,-2] == '.el' ]]; then
		export LD_LIBRARY_PATH=$HOME/sw/lib:$LD_LIBRARY_PATH
	else

	fi
	export PYTHONPATH="$PYTHONPATH:$HOME/dropbox/pymol"
	export PYTHONPATH=$HOME/Dropbox/lib/python:$PYTHONPATH

	#export ROSETTA3_DB=$HOME/rosetta/dev/database
	export EDITOR='subl -w'

__git_files () {
	_wanted files expl 'local files' _files
}

function current_branch {
	# git branch | grep \* | cut -b 3-
}

################################################# UTILS ################################################################
	function biounit {
		open -a /Applications/MacPyMOL.app /data/pdb/biounit/${1:2:2}/${1:1:4}.pdb*.gz
	}
	function pdbxml  { gzcat /data/pdb/XML-noatom/${1:2:2}/${1:1:4}-noatom.xml.gz > ~/tmp/$1.xml; subl ~/tmp/$1.xml };
	function pdbbbonly {
		for f in $@; do
		echo removing not N+CA+C+O+CB $f
		egrep '  N   |  CA  |  C   |  CB  |  O   ' $f > .tmp_pdbbbonly
		mv .tmp_pdbbbonly $f
		done
	}
	function hackfix2b98 {
		for i in $@; do
			echo $i
			egrep -v "^ATOM      .  ... THR A   1   " $i > $i.tmp
			mv $i.tmp $i
		done
	}

################################################ ALIASES ################################################################
	alias xgr='xargs grep'
	alias fxgr='find . -type f | xargs grep'
	alias fxgrhh='find . -type f -name \*.hh | xargs grep'
	alias fxgrcc='find . -type f -name \*.cc | xargs grep'
	alias rso="cd $HOME/rosetta/dev/source/cmake/build_omp"
	alias rsd="cd $HOME/rosetta/dev/source/cmake/build_debug"
	alias rss="cd $HOME/rosetta/dev/source"
	alias cmakeupmy="pushd ../..; ./update_options.sh; pushd cmake; ./make_project.py my; popd; popd; ninja -k0 -j8"
	alias cmakeupall="pushd ../..; ./update_options.sh; pushd cmake; ./make_project.py all; popd; popd; ninja -k0 -j8"
	alias clangomp="clang++ -I/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/include -ferror-limit=1"
	alias extractchainAB=" egrep '^...... .... .... ... (A|B) ... '"


source $HOME/.vars
source $HOME/.zshrc.this_machine
