# Optimizing Square Matrix Symmetry Verification and Transposition

Optimizing Square Matrix Symmetry Verification and Transposition Using Implicit Parallelization Techniques, OpenMP and MPI

## Assignment links
For H1 (H1 Folder) :  [H1 (Implicit Parallelization and OpenMP)](#h1---implicit-parallelization-and-openmp) <br>
    - [Requirements](#requirements) <br>
    - [Project Structure](#project-structure) <br>
    - [Instructions for Reproduction](#instruction-for-reproduction) <br>
[You can find the H1 Paper here](H1/H1.pdf) <br>
    

For H2 (H2 Folder) :  [H2 (OpenMP and MPI)](#h2---mpi-and-openmp) <br>
    - [Requirements](#requirements-1) <br>
    - [Project Structure](#project-structure-1) <br>
    - [Instructions for Reproduction](#instruction-for-reproduction-1) <br>
[You can find the H2 Paper here](H2/H2.pdf) <br>

***

# H1 - Implicit Parallelization and OpenMP

The folders contains the c codes and some python scripts to simplify and automatize the data elaboration and procurement.
Each folder represents a version: Imp* folders represent optimized sequential versions, Omp* folders represent parallelized versions with OpenMP.
Inside each version there is a "Data" folder containing all related data: *_processed.csv contains the data for Transposition and Check symmetry;  *_processed_SpeeEff.csv, if present, collects data related to Speedup and Efficiency (only for omp versions).
Inside each Omp* folder there is a "Charts" folder with graphs related to that version.

## Requirements

- Gcc 9.1
- UniTn HPC
- OpenMP 3.1
- Python 3.7.2

## Project Structure

```
H1/
├── 64vs32vs16/ -> finding the best block size
├── ImpAll/ -> sequential optimized with the best configuration
├── ImpOnlyO2/ -> sequential optimized with only O2 flag
├── ImpOnlyFlags/ -> sequential optimized with only compilation flags without O2
├── ImpPragma/ -> sequential optimized with only pragma instructions
├── OmpDynamic/ -> parallel with Blocking and schedule(Dynamic) - Best configuration
├── OmpNoBlocking/ -> parallel without Blocking
├── OmpStatic/ -> parallel with Blocking and schedule(Static)
├── BW/ -> contains all data about the Bandwidth
├── H1.pdf -> a pdf version of the paper of the project
└── scriptTotal.pbs
```

## Instructions for Reproduction 

To run the script:
```bash
cd H1
```
Here edit line 17 of the file scriptTotal.pbs: cd /home/name.surname/H1/64vs32vs16 substituting name.surname with the credential of your account and then:
```bash
qsub scriptTotal.pbs
```
It requires some time to complete (only the run takes around 40 minutes).

In each folder, a *_processed.csv file will be created with all data
and, if the folder contains an OpenMP version, a new *_processed_SpeeEff.csv file will be created containing Speedup and Efficiency data

-----

Graphs available only for OpenMP versions

To generate graphs, firtly load the python module:
```bash
module load python-3.7.2
```

Then enter the desired version folder and run:
```python
python3 plot.py <*_processed.csv>
```
For logarithmic scale graphs:
```python
python3 plot_log.py <*_processed.csv>
```
For speedup and efficiency graphs:
```python
python3 plot_eff.py <*_processed_SpeeEff.csv>
```
For logarithmic scale speedup and efficiency graphs:
```python
python3 plot_eff_log.py <*_processed_SpeeEff.csv>
```
To calculate minimum times for Bandwidth calculation, enter the OmpDynamic folder and run:
```python
python3 Find_Min.py <*_processed.csv>
```
This will create a Min_MatTraspose.csv file with times by Dim.
For bandwidth calculation, formulas must then be applied manually.

---

---

# H2 - MPI and OpenMP

The folder contain the c codes and some python scripts to simplify and automatize the data elaboration and procurement.
Inside each folder there is a "Data" folder containing all related data: *_processed.csv contains the data for Transposition and Check symmetry; *_processed_SpeeEff.csv, if present, collects data related to Speedup, Efficiency and Scaling.
In some folders there is: a "Charts" folder with all the graphs, a "Scaling" and a "Speedup&Efficiency" folder that contain the performance data and graphs of it
In mpiBlock folder we can find the algorithms that divide the matrix in block, sent the blocks among processes that transpose it and sends back to the root process

## Requirements

- Gcc 9.1
- UniTn HPC node 07 (or another 96 cores)
- Mpich 3.2.1
- OpenMP 3.1
- Python 3.7.2

## Project Structure
```
H2/
├── BaseVsAlltoallVsDatatype/ -> confront among different MPI versions
|       ├── MPI_ALL.c -> represent an MPI version using all to all
|       ├── MPI_Base.c -> represent an MPI version using only send and recieve
|       └── MPI_DataT.c -> represent an MPI version using MPI Datatypes
├── comp/ -> comparison between MPI and OMP and SEQ best versions
|       ├── OMP
|       ├── MPI 
|       └── SEQ
├── mpiBlock/ -> comparison between the block dividing MPI version (parallel and sequential) vs normal MPI version
|       └── BlockComp/ -> comparison between different MPI block dividing versions to find the best
|              ├── MPIB_Basic.c -> basic version 
|              ├── MPIB_DataTV.c -> using MPI Datatypes and Variable communicators
|              ├── MPIB_nvnd.c -> only using
|              └── MPIB_onlyD.c -> only using datatypes
├── info.txt -> the UniTN HPC node info
├── H2.pdf -> a pdf version of the paper of the project
└── script.pbs
```

## Instructions for Reproduction 

To run the script:
```bash
cd H2
```
Here edit line 26 of the file script.pbs: cd /home/name.surname/H2 substituting name.surname with the credential of your account, if you want to receive a notification via email when the simulation starts/finish change, in the line 13, "mose.arcaro@studenti.unitn.it" with your email and then:
```bash
qsub script.pbs
```
It requires some time to complete (only the run takes around 1 hour and 50 minutes).

In each folder, a *_processed.csv file will be created with all data 
and, if is a comparison folder, a new *_processed_SpeeEff.csv file will be created containing Speedup and Efficiency data.
In Scaling folder, if present, a new *_processed_with_scaling.csv will be created containing Scaling data

-----

Graphs available only for "comp" and "mpiBlock" folders

To generate graphs, firtly load the python module:
```bash
module load python-3.7.2
```

Then enter the desired version folder 
```bash
cd <desired folder>
```
and run:
```python
python3 plot.py <*_processed.csv>
```
For logarithmic scale graphs:
```python
python3 plot_log.py <*_processed.csv>
```
For speedup and efficiency graphs:
```python
python3 plot_eff.py <*_processed_SpeeEff.csv>
```
For logarithmic scale speedup and efficiency graphs:
```python
python3 plot_eff_log.py <*_processed_SpeeEff.csv>
```
For scaling graphs enter the designed directory:
```bash
cd Scaling
```
and run the script:
```python
python3 plotScal.py <*_processes_with_scaling.csv>
```

now a new "grafici_scaling" floder will be created with all scaling graphs




