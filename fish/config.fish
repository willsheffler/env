


eval /home/sheffler/miniconda3/bin/conda "shell.fish" "hook" $argv | source
conda activate rfd
# source $HOME/venv/sci/bin/activate.fish

set -x __NV_PRIME_RENDER_OFFLOAD 1
set -x __GLX_VENDOR_LIBRARY_NAME nvidia

set -x OMP_NUM_THREADS 4
set -x MKL_NUM_THREADS 4

#set -x LD_PRELOAD '/home/sheffler/tools/graphbolt/lib/mimalloc/out/release/libmimalloc.so'
set -x APPTAINER_CONTAINER '/software/containers/users/sheffler/rf_diffusion_aa_py310.sif'
set -x PY /home/sheffler/venv/sci/bin/python
# function fish_prompt; echo "-> "; end
# function fish_prompt_right; echo ""; end
