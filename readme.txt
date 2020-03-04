****main commands to run****
nvcc -o op main2.cu
./op <fullgraph file name> <SSSP file name> <changeEdges file name> <no. of nodes> <no. of edges>





****How to run****
1. command to access gpu node:
srun --partition=gpu --gres=gpu --mem=4gb --ntasks-per-node=2 --nodes=1 --pty

2. command to compile 
nvcc -o op main2.cu

3.create the job batch.sub 
#!/bin/sh
#SBATCH --time=03:15:00
#SBATCH --mem-per-cpu=1024
#SBATCH --job-name=cuda
#SBATCH --partition=gpu
#SBATCH --gres=gpu
#SBATCH --error=/work/[groupname]/[username]/job.%J.err
#SBATCH --output=/work/[groupname]/[username]/job.%J.out

module load cuda
./op.exe


4. run the job 
sbatch batch.sub


****Debug using gdb****
1. gdb op
2. run
3. bt
