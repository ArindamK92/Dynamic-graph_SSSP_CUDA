****Input Files****
You need to supply changeEdges.txt, fullGraph.txt and SSSP.txt files as input. The files given here are samples. Real files will contain the actual inputs folling the same format.
keep changeEdges.txt, fullGraph.txt and SSSP.txt in the same folder where you will keep main1.cu


****How to run****
1. command to access gpu node:
srun --partition=gpu --gres=gpu --mem=4gb --ntasks-per-node=2 --nodes=1 --pty

2. command to compile 
nvcc -o op main1.cu

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