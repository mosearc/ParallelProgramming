# Optimizing Square Matrix Symmetry Verification and Transposition

Optimizing Square Matrix Symmetry Verification and Transposition Using Implicit Parallelization Techniques and OpenMP  
Ogni cartella rappresenta una versione e dentro ognuna di essa è presente la cartella "Data" con i dati e grafici di quella versione
(grafici disponibili solo per le versioni con OpenMP
## Requisiti

- Gcc 9.1
- UniTn HPC
- OpenMP 3.1
- Pyhton 3.7.2

## Struttura del Progetto

```
Total/
├── ImpAll/ -> sequential Optimized with the best configuration
├── ImpOnliO2/ -> sequential optimized with only O2 flag
├── ImpOnlyFlags/ -> sequantial optimized with only compilation flags without O2
├── ImpPragma/ -> sequential optimized with only pragma instructions
├── OmpDynamic/ -> parallel with Blocking and schedule(Dynamic)
├── OmpNoBlocking/ -> parallel without Blocking
├── OmpStatic/ -> parallel with Blocking and schedule(Static)
├── BW.csv/ -> Data for Bandwidth
└── scriptTotal.pbs
```


## Avvio

Per avviare in modalità sviluppo:
```bash
qsub scriptTotal.pbs
```
in ogni cartella verrà creato un file .._processed.csv con tutti i dati
e, se la cartella contine una versione OpenMP verrà creato un nuovo file .._processed_SpeeEff.csv che contiene anche i dati degli speedup ed efficiency

se si vuole generare i grafici entrare nella cartella della versione desiderata ed eseguire:
```python
python3 plot.py <file_processed.csv>
```
se si vuole i grafici in scala logaritmica entrare nella cartella della versione desiderata ed eseguire:
```python
python3 plot_log.py <file_processed.csv>
```
se si vogliono i grafici di speedup ed efficiency entrare nella cartella della versione desiderata ed eseguire:
```python
python3 plot_eff.py <file_processed_SpeeEff.csv>
```

se si vogliono i grafici di speedup ed efficiency in scala logaritmica entrare nella cartella della versione desiderata ed eseguire:
```python
python3 plot_eff_log.py <file_processed_SpeeEff.csv>
```
se si vogliono calcolare i tempi minimi per calcolare la bandwidth entrare nella cartella OmpDynamic ed eseguire:
```python
python3 Find_Min.py <file_processed*.csv>
```
per il calcolo della Bandwidth è necessario poi appliare le formule a mano

