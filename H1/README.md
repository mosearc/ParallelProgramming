# Optimizing Square Matrix Symmetry Verification and Transposition

Optimizing Square Matrix Symmetry Verification and Transposition Using Implicit Parallelization Techniques and OpenMP 

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
Total/
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

## Execution

To run the script:
```bash
cd H1
```
Here edit line 17 of the file scriptTotal.pbs: cd /home/name.surname/Total/64vs32vs16 substituting name.surname with the credential of your account and then:
```bash
qsub scriptTotal.pbs
```
It requires some time to complete (around 40 minutes).

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
