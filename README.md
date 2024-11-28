# Optimizing Square Matrix Symmetry Verification and Transposition

Optimizing Square Matrix Symmetry Verification and Transposition Using Implicit Parallelization Techniques and OpenMP  
Ogni cartella contiene il proprio readme con la descrizione del contenuto

## Requisiti

- Gcc 9.1
- UniTn HPC
- OpenMP 3.1
- Pyhton 3.7.2

## Avvio

Per avviare in modalità sviluppo:
```bash
qsub script.pbs
```
in ogni cartella verrà creato un file .._processed.csv con tutti i dati
e, se la cartella contine una versione OpenMP verrà creato un nuovo file .._processed_SpeeEff.csv che contiene anche i dati degli speedup ed efficiency

se si vuole generare i grafici eseguire:
```python
python3 plot.py <file_processed.csv>
```
se si vuole i grafici in scala logaritmica eseguire:
```python
python3 plot_log.py <file_processed.csv>
```
se si vogliono i grafici di speedup ed efficiency eseguire:
```python
python3 plot_eff.py <file_processed_SpeeEff.csv>
```

se si vogliono i grafici di speedup ed efficiency in scala logaritmica eseguire:
```python
python3 plot_eff_log.py <file_processed_SpeeEff.csv>
```
se si vogliono calcolare i tempi minimi per calcolare la bandwidth eseguire:
```python
python3 Find_Min.py <file_processed*.csv>
```
in output si avra un file con le dimenisioni e la corrispettiva tempistica, eseguire poi il calcolo a mano per trovare la bandwidth


