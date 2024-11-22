#!/bin/bash
#SBATCH --job-name=KATZ
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=gpu

# Run the solution experiment on a single CPU/GPU node.
# The results are redirected to logfiles in the katzsolution folder.
# - There is one log file for each matrix and for each value of alpha.
# - Objective of this test is showing that there is a reliable procedure for
#   performing matrix-vector products with the modified Laplacian (Katz)


# Load the module for the Toeplitz machine:
module purge
module load gpu-gcc/12.2.0 gpu-metis/5.1.0-gcc-12.2.0 gpu-cuda/12.3.1-gcc-12.2.0 gpu-openmpi/4.1.6-cuda-12.3.1-gcc-12.2.0 gpu-openblas/0.3.26-gcc-12.2.0
module list

# Set the environment variable for the number of OPENMP processes (needs to be equal to --cpus-per-task)
export OMP_NUM_THREADS=64

matrixlist=('../testmatrices/mm/PGPgiantcompo.mtx' '../testmatrices/mm/USpowerGrid.mtx' '../testmatrices/mm/belgium_osm.mtx' '../testmatrices/mm/ca-CondMat.mtx' '../testmatrices/mm/ca-HepPh.mtx' '../testmatrices/mm/dictionary28.mtx' '../testmatrices/mm/human_gene1.mtx' '../testmatrices/mm/human_gene2.mtx' '../testmatrices/mm/roadNet-CA.mtx' '../testmatrices/mm/roadNet-PA.mtx' '../testmatrices/mm/roadNet-TX.mtx' '../testmatrices/mm/usroads-48.mtx')
declare -A alphavalue
alphavalue[0,0]=0.00002437
alphavalue[0,1]=0.00016605
alphavalue[0,2]=0.00113131
alphavalue[0,3]=0.00770750
alphavalue[1,0]=0.00016061
alphavalue[1,1]=0.00109421
alphavalue[1,2]=0.00745475
alphavalue[1,3]=0.05078861
alphavalue[2,0]=0.00038737
alphavalue[2,1]=0.00263912
alphavalue[2,2]=0.01798009
alphavalue[2,3]=0.12249694
alphavalue[3,0]=0.00002789
alphavalue[3,1]=0.00019001
alphavalue[3,2]=0.00129449
alphavalue[3,3]=0.00881929
alphavalue[4,0]=0.00000410
alphavalue[4,1]=0.00002795
alphavalue[4,2]=0.00019042
alphavalue[4,3]=0.00129733
alphavalue[5,0]=0.00006224
alphavalue[5,1]=0.00042405
alphavalue[5,2]=0.00288901
alphavalue[5,3]=0.01968256
alphavalue[6,0]=0.00000030
alphavalue[6,1]=0.00000204
alphavalue[6,2]=0.00001391
alphavalue[6,3]=0.00009480
alphavalue[7,0]=0.00000033
alphavalue[7,1]=0.00000227
alphavalue[7,2]=0.00001547
alphavalue[7,3]=0.00010542
alphavalue[8,0]=0.00030109
alphavalue[8,1]=0.00205131
alphavalue[8,2]=0.01397540
alphavalue[8,3]=0.09521328
alphavalue[9,0]=0.00032156
alphavalue[9,1]=0.00219080
alphavalue[9,2]=0.01492572
alphavalue[9,3]=0.10168775
alphavalue[10,0]=0.00028091
alphavalue[10,1]=0.00191382
alphavalue[10,2]=0.01303869
alphavalue[10,3]=0.08883158
alphavalue[11,0]=0.00036175
alphavalue[11,1]=0.00246460
alphavalue[11,2]=0.01679112
alphavalue[11,3]=0.11439657


i=0
for matrix in ${matrixlist[@]}
do
  for j in 0 1 2 3
  do
  alpha=${alphavalue[${i},${j}]}
  echo "Working on matrix ${matrix} for α = ${alpha}"

srun ./amg_dkatz 2>&1 > katzsolution/log_cpu_matrix${i}_alpha${j}.txt << EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
${matrix}                   ! Other matrices from: http://math.nist.gov/MatrixMarket/ or
NONE                        ! rhs          ! http://www.cise.ufl.edu/research/sparse/matrices/index.html
NONE                        ! Initial guess
NONE                        ! Reference solution
MM                          ! File format: MatrixMarket or Harwell-Boeing
CSR                         ! Storage format: CSR COO JAD
GRAPH			                  ! PART (partition method): BLOCK GRAPH
${alpha}                    ! α value
1.0                         ! μ value
FCG                         ! Iterative method: BiCGSTAB BiCGSTABL BiCG CG CGS FCG GCR RGMRES
2                           ! ISTOPC
00100                       ! ITMAX
1                           ! ITRACE
30                          ! IRST (restart for RGMRES and BiCGSTABL)
1.d-9                       ! EPS
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-SVCYC-L1JAC1-SOC1-L1JAC  ! Longer descriptive name for preconditioner (up to 20 chars)
ML                          ! Preconditioner type: NONE JACOBI GS FBGS BJAC AS ML
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Smoother type JACOBI FBGS GS BWGS BJAC AS. For 1-level, repeats previous.
1                           ! Number of sweeps for smoother
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
ILU                         ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
0                           ! Fill level P for ILU(P) and ILU(T,P)
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                        ! Second (post) smoother, ignored if NONE
1                           ! Number of sweeps for (post) smoother
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
ILU                         ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
0                           ! Fill level P for ILU(P) and ILU(T,P)
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                      ! Type of multilevel CYCLE: VCYCLE WCYCLE KCYCLE MULT ADD
1                           ! Number of outer sweeps for ML
6                           ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                          ! Target coarse matrix size; if <0, lib default
SMOOTHED                    ! Type of aggregation: SMOOTHED UNSMOOTHED
DEC                         ! Parallel aggregation: DEC, SYMDEC
SOC1                        ! Algorithm
8                           ! Size of the aggregates
NATURAL                     ! Ordering of aggregation NATURAL DEGREE
-1.5                        ! Coarsening ratio, if < 0 use library default
FILTER                      ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0                   ! Smoothed aggregation threshold, ignored if < 0
-2                          ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025                  ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
UMF                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
DIST                        ! Coarsest-level matrix distribution: DIST  REPL
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
EOF

srun ./amg_dkatz 2>&1 > katzsolution/log_gpu_matrix${i}_alpha${j}.txt << EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
${matrix}                   ! Other matrices from: http://math.nist.gov/MatrixMarket/ or
NONE                        ! rhs          ! http://www.cise.ufl.edu/research/sparse/matrices/index.html
NONE                        ! Initial guess
NONE                        ! Reference solution
MM                          ! File format: MatrixMarket or Harwell-Boeing
CSRG                        ! Storage format: CSR COO JAD
GRAPH			                  ! PART (partition method): BLOCK GRAPH
${alpha}                    ! α value
1.0                         ! μ value
FCG                         ! Iterative method: BiCGSTAB BiCGSTABL BiCG CG CGS FCG GCR RGMRES
2                           ! ISTOPC
00100                       ! ITMAX
1                           ! ITRACE
30                          ! IRST (restart for RGMRES and BiCGSTABL)
1.d-9                       ! EPS
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-SVCYC-L1JAC1-SOC1-L1JAC  ! Longer descriptive name for preconditioner (up to 20 chars)
ML                          ! Preconditioner type: NONE JACOBI GS FBGS BJAC AS ML
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Smoother type JACOBI FBGS GS BWGS BJAC AS. For 1-level, repeats previous.
1                           ! Number of sweeps for smoother
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
ILU                         ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
0                           ! Fill level P for ILU(P) and ILU(T,P)
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                        ! Second (post) smoother, ignored if NONE
1                           ! Number of sweeps for (post) smoother
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
ILU                         ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
0                           ! Fill level P for ILU(P) and ILU(T,P)
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                      ! Type of multilevel CYCLE: VCYCLE WCYCLE KCYCLE MULT ADD
1                           ! Number of outer sweeps for ML
6                           ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                          ! Target coarse matrix size; if <0, lib default
SMOOTHED                    ! Type of aggregation: SMOOTHED UNSMOOTHED
DEC                         ! Parallel aggregation: DEC, SYMDEC
SOC1                        ! Algorithm
8                           ! Size of the aggregates
NATURAL                     ! Ordering of aggregation NATURAL DEGREE
-1.5                        ! Coarsening ratio, if < 0 use library default
FILTER                      ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0                   ! Smoothed aggregation threshold, ignored if < 0
-2                          ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025                  ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
UMF                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
DIST                        ! Coarsest-level matrix distribution: DIST  REPL
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
EOF

done
let i=i+1
done
