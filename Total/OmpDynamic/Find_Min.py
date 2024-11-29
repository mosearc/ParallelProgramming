import pandas as pd
import sys
import os

def find_min_rolling_mean_mattrans(input_file):
    # Load the dataset
    data = pd.read_csv(input_file)

    # Filtra i dati per 'Version' == 'OMP'
    omp_data = data[data['Version'] == 'OMP']

    # Raggruppa per 'Dim' e calcola il minimo di 'RollingMean_MatTrans' per ogni gruppo
    min_values = omp_data.groupby('Dim')['RollingMean_MatTranspose'].min()

    # Stampa i risultati
    print("Minimi valori di RollingMean_MatTranspose per ogni Dim (Version = OMP):")
    print(min_values)

    # Se vuoi salvare i risultati in un file CSV
    output_file = 'Min_MatTranspose.csv'
    min_values.to_csv(output_file)
    print(f"Risultati salvati in {output_file}")

if __name__ == "__main__":
    # Verifica se lo script è stato chiamato con l'argomento giusto
    if len(sys.argv) != 2:
        print("Uso: python script.py <file_csv_input>")
        sys.exit(1)

    # Ottieni il nome del file da linea di comando
    input_file = sys.argv[1]

    # Verifica se il file esiste
    if not os.path.isfile(input_file):
        print(f"Errore: Il file {input_file} non esiste.")
        sys.exit(1)

    # Trova i minimi di RollingMean_MatTrans per ogni Dim
    find_min_rolling_mean_mattrans(input_file)

