#!/bin/bash

#SBATCH --job-name=phantom_pcg
#SBATCH --output=phantom_pcg_%j.out
#SBATCH --error=phantom_pcg_%j.err
#SBATCH --nodelist=bigMem0
#SBATCH --time=18:00:00
#SBATCH --exclusive

export PATH=$PATH:/home/jyliu/NUFFT2D/extern/fftw-3.3.10/build_double/bin:/home/jyliu/NUFFT2D/extern/fftw-3.3.10/build_single/bin
export LIBRARY_PATH=$LIBRARY_PATH:/home/jyliu/NUFFT2D/extern/fftw-3.3.10/build_double/lib:/home/jyliu/NUFFT2D/extern/fftw-3.3.10/build_single/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jyliu/NUFFT2D/extern/fftw-3.3.10/build_double/lib:/home/jyliu/NUFFT2D/extern/fftw-3.3.10/build_single/lib
export CPATH=$CPATH:/home/jyliu/NUFFT2D/extern/fftw-3.3.10/build_double/include:/home/jyliu/NUFFT2D/extern/fftw-3.3.10/build_single/include

module unload MATLAB
module load MATLAB/R2023b
matlab -r 'cd /home/jyliu/NUFFT2D; nufft2d_startup; cd experiments/typeII_2d_phantom_pcg; exp_typeII_2d_phantom_pcg;'
