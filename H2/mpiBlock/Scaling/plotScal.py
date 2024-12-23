import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import numpy as np

def plot_and_save_scaling(input_file, output_dir):
    # Caricare il file CSV
    df = pd.read_csv(input_file)

    # Controllare che le colonne necessarie esistano
    required_columns = ['Processes', 'scaling_RollingMean_MatTranspose']
    if not all(col in df.columns for col in required_columns):
        raise ValueError("Il file deve contenere le colonne 'Processes' e 'scaling_RollingMean_MatTranspose'.")

    # Creare la cartella di output se non esiste
    os.makedirs(output_dir, exist_ok=True)

    # Identificare i punti di partenza (Processes == 1)
    start_indices = df[df['Processes'] == 1].index
    start_indices = list(start_indices) + [df.index[-1] + 1]

    # Generare grafici per ogni segmento con titolo basato su potenze di 2
    current_power = 16

    for i in range(len(start_indices) - 1):
        start, end = start_indices[i], start_indices[i+1]
        segment = df.iloc[start:end]

        # Creare una linea di decremento inverso (1 / x)
        x_values = segment['Processes']
        inverse_bg = [1 / x if x != 0 else None for x in x_values]

        # Grafico con `scaling_RollingMean_MatTranspose` e decremento inverso
        plt.figure(figsize=(12, 7))
        plt.plot(segment['Processes'], segment['scaling_RollingMean_MatTranspose'], marker='o', label='scaling_RollingMean_MatTranspose')
        plt.plot(segment['Processes'], inverse_bg, linestyle='--', color='gray', label='Decremento Inverso (1/x)')

        plt.title(f'Scaling RollingMean_MatTranspose - Potenza di 2: {current_power}')
        plt.xlabel('Processes')
        plt.ylabel('Scaling Values')
        plt.legend()
        plt.grid(True)

        # Salvare il grafico
        output_path = os.path.join(output_dir, f'scaling_potenza_{current_power}.png')
        plt.savefig(output_path)
        plt.close()

        print(f"Grafico salvato: {output_path}")

        # Incrementare la potenza di 2
        current_power *= 2


if __name__ == "__main__":
    # Parser degli argomenti da riga di comando
    parser = argparse.ArgumentParser(description='Genera e salva grafici per "scaling_RollingMean_MatTranspose" con decremento inverso (1/x).')
    parser.add_argument('input_file', type=str, help='Percorso del file CSV di input')
    parser.add_argument('--output_dir', type=str, default='scaling_graphs', help='Cartella di output per i grafici (default: scaling_graphs)')

    args = parser.parse_args()

    # Esegui la funzione principale
    plot_and_save_scaling(args.input_file, args.output_dir)
