# Optimizing Square Matrix Symmetry Verification and Transposition

Optimizing Square Matrix Symmetry Verification and Transposition Using Implicit Parallelization Techniques and OpenMP 


Ogni cartella rappresenta una versione: le cartelle Imp* rappresentano le versioni sequenziali ottimizzate, le cartelle Omp* rappresentano le versioni parallelizzate con OpenMP<br />
Dentro ogni versione è presente la cartella "Data" dove è possibile trovare tutti i relativi dati: *_processed_SpeeEff.csv raccoglie i dati relativi agli Speedup ed Efficency (solo per versioni omp)<br />
Dentro ogni cartella Omp* è presente la cartella "Charts" con i grafici relativi a quella versione.

## Requisiti

- Gcc 9.1
- UniTn HPC
- OpenMP 3.1
- Pyhton 3.7.2

## Struttura del Progetto

```
Total/
├── ImpAll/ -> sequential optimized with the best configuration
├── ImpOnlyO2/ -> sequential optimized with only O2 flag
├── ImpOnlyFlags/ -> sequantial optimized with only compilation flags without O2
├── ImpPragma/ -> sequential optimized with only pragma instructions
├── OmpDynamic/ -> parallel with Blocking and schedule(Dynamic) - Best configurtion
├── OmpNoBlocking/ -> parallel without Blocking
├── OmpStatic/ -> parallel with Blocking and schedule(Static)
├── BW.csv/ -> Data for Bandwidth
└── scriptTotal.pbs
```


## Avvio

Per avviare lo script (it requires around 40 minutes to complete):
```bash
cd Total
qsub scriptTotal.pbs
```
in ogni cartella verrà creato un file *_processed.csv con tutti i dati
e, se la cartella contine una versione OpenMP, verrà creato un nuovo file *_processed_SpeeEff.csv che contiene anche i dati degli Speedup ed Efficiency



Grafici disponibili solo per versioni OpenMP

se si vuole generare i grafici entrare nella cartella della versione desiderata ed eseguire:
```python
python3 plot.py <*_processed.csv>
```
se si vuole i grafici in scala logaritmica entrare nella cartella della versione desiderata ed eseguire:
```python
python3 plot_log.py <*_processed.csv>
```
se si vogliono i grafici di speedup ed efficiency entrare nella cartella della versione desiderata ed eseguire:
```python
python3 plot_eff.py <*_processed_SpeeEff.csv>
```
se si vogliono i grafici di speedup ed efficiency in scala logaritmica entrare nella cartella della versione desiderata ed eseguire:
```python
python3 plot_eff_log.py <*_processed_SpeeEff.csv>
```
se si vogliono calcolare i tempi minimi per calcolare la Bandwidth entrare nella cartella OmpDynamic ed eseguire, verrà creato un file Min_MatTraspose.csv con i tempi per Dim:
```python
python3 Find_Min.py <*_processed.csv>
```
per il calcolo della Bandwidth è necessario poi appliare le formule a mano

