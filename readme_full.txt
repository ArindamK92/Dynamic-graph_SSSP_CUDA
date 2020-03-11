How to create changeedge file:
____________________________________
nvcc -o op_createChangeEdges createChangedEdges.cpp
OR we can use g++ -o createCE createChangedEdges.cpp -std=c++11 -O3

op_createChangeEdges.exe <fullGraph.txt> <no. of nodes> <no of change edge> <percentage of insertion>
OR
./createCE <fullGraph.txt> <no. of nodes> <no of change edge> <percentage of insertion>

How to create SSSP input file(sequential):
___________________________________________
nvcc -o output <SeqSSSPmain.cpp>
or
g++ -o output <filename mainSSSP.cpp> -std=c++11 -O3

./output <fullgraphfilename> <no. of nodes>


How to connect disconnected nodes:
_____________________________________________________
nvcc -o op_createChangeEdges.exe createChangedEdges.cpp
op_connectDisconnectedNodes.exe <SSSPfile> > <newFilename to store new edges>

add the new synthesized edges in full graph and again find seqSSSP. Use new updated fullgraph and new seqSSSP file for dynamic SSSP.


How to run parallel CUDA SSSP code:
____________________________________
****main commands to run****
nvcc -o op main2.cu
./op <fullgraph file name> <SSSP file name> <changeEdges file name> <no. of nodes> <no. of edges * 2 (or total number of lines in fullgraph file)>


Forge cuda commands:
sinteractive -p cuda --time=01:00:00 --gres=gpu:1 --mem=24000M
module load cuda/10.1 gnu/7.2.0
nvcc -o output main1.cu -std=c++11 -O3
./output


****How to run(for Sriram's cluster)****
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

To get total lines in a file in linux:
wc -l fullGraph.txt

How to get nvidia gpu details:
nvidia-smi
