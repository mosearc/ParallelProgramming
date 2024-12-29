import pandas as pd
import argparse
import os

def calculate_scaling(input_file, output_file):
    # Carica il file CSV
    df = pd.read_csv(input_file)

    # Controlla che le colonne necessarie esistano
    required_columns = ['Processes', 'RollingMean_MatTranspose']
    if not all(col in df.columns for col in required_columns):
        raise ValueError("Il file deve contenere le colonne 'Processes' e 'RollingMean_MatTranspose'.")

    # Calcola la colonna "scaling"
    scaling = []
    scaling_value = None

    for index, row in df.iterrows():
        if row['Processes'] == 1:
            # Imposta il nuovo valore di riferimento per il calcolo dello scaling
            scaling_value = row['RollingMean_MatTranspose']
        if scaling_value is not None:
            scaling.append(scaling_value / row['RollingMean_MatTranspose'])
        else:
            scaling.append(None)

    df['scaling_RollingMean_MatTranspose'] = scaling

    # Salva il risultato in un nuovo file
    df.to_csv(output_file, index=False)
    print(f"File salvato con successo come: {output_file}")


if __name__ == "__main__":
    # Parser degli argomenti da riga di comando
    parser = argparse.ArgumentParser(description='Calcola la colonna "scaling" in un file CSV.')
    parser.add_argument('input_file', type=str, help='Percorso del file CSV di input')
    parser.add_argument('--output', type=str, default=None, help='Percorso del file CSV di output (opzionale)')

    args = parser.parse_args()

    # Definisci il file di output se non è specificato
    output_file = args.output if args.output else os.path.splitext(args.input_file)[0] + '_with_scaling.csv'

    # Esegui la funzione principale
    calculate_scaling(args.input_file, output_file)
