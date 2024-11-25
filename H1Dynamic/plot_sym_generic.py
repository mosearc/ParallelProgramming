import argparse
import pandas as pd
import matplotlib.pyplot as plt


#python3 plot_sym_generic.py outputMeansLast.csv



# Configurazione di argparse per accettare il file da riga di comando
parser = argparse.ArgumentParser(description="Crea grafici Sym in base ai dati forniti.")
parser.add_argument("file", type=str, help="Il file CSV contenente i dati")
args = parser.parse_args()

# Caricamento del file CSV
try:
    data = pd.read_csv(args.file)
except FileNotFoundError:
    print(f"Errore: Il file '{args.file}' non esiste.")
    exit(1)
except pd.errors.EmptyDataError:
    print(f"Errore: Il file '{args.file}' è vuoto o non è valido.")
    exit(1)

# Assicurati che il file abbia le colonne richieste
required_columns = {'Dim', 'Version', 'Threads', 'Sym'}
if not required_columns.issubset(data.columns):
    print(f"Errore: Il file deve contenere le colonne: {', '.join(required_columns)}")
    exit(1)

# Creazione dei grafici SYM
unique_dims = data['Dim'].unique()
for dim in unique_dims:
    subset = data[data['Dim'] == dim]
    plt.figure(figsize=(8, 6))
    
    # Traccia per la versione OMP (colore distintivo)
    omp_data = subset[subset['Version'] == 'OMP']
    if not omp_data.empty:
        plt.plot(omp_data['Threads'], omp_data['Sym'], marker='o', color='blue', label='Open MP')
    
    # Linea costante per la versione SEQ (colore distintivo)
    seq_data = subset[subset['Version'] == 'SEQ']
    if not seq_data.empty and not seq_data['Sym'].isna().all():
        mean_seq_sym = seq_data['Sym'].mean()
        plt.axhline(y=mean_seq_sym, linestyle='--', color='green', label='Sequential')
    
    # Linea costante per la versione IMP (colore distintivo)
    imp_data = subset[subset['Version'] == 'IMP']
    if not imp_data.empty and not imp_data['Sym'].isna().all():
        mean_imp_sym = imp_data['Sym'].mean()
        plt.axhline(y=mean_imp_sym, linestyle='--', color='orange', label='Implicit')
    
    # Impostazioni del grafico
    plt.title(f'Check simmetry vs Threads per Dim ({dim})')
    plt.xlabel('Threads [#]')
    plt.ylabel('Time [s]')
    plt.legend()
    plt.grid(True)
    
        # Salvataggio del grafico
    filename = f"Sym_vs_Threads_Dim_{dim}.png"
    plt.savefig(filename, format='png')  # Salva il grafico come file PNG
    print(f"Grafico salvato come: {filename}")
    
    plt.show()

