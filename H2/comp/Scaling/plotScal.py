import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def plot_scaling_trends(file_path):
    # Caricare il file CSV
    df = pd.read_csv(file_path)

    # Creare una cartella per salvare i grafici
    output_dir = 'grafici_scaling_dual'
    os.makedirs(output_dir, exist_ok=True)

    # Aggiungere colonna per rilevare i cambi di versione
    df['Version_Change'] = df['Version'] != df['Version'].shift()
    df['Group'] = df['Version_Change'].cumsum()

    # Raggruppare per gruppi distinti di versioni
    version_groups = df.groupby('Group')
    grouped_versions = list(version_groups)

    # Creare grafici a coppie con entrambe le metriche e la linea 1/x
    for i in range(0, len(grouped_versions), 2):
        plt.figure(figsize=(12, 7))

        for j in range(2):
            if i + j < len(grouped_versions):
                group_id, group = grouped_versions[i + j]
                version = group['Version'].iloc[0]

                # Tracciare entrambe le metriche
                plt.plot(group['Processes'], group['scaling_RollingMean_MatTranspose'],
                         marker='o', label=f'{version} - MatTranspose ')
                plt.plot(group['Processes'], group['scaling_RollingMean_CheckSym'],
                         marker='s', label=f'{version} - CheckSym ')

        # Aggiungere la linea di riferimento 1/x
        x_values = np.linspace(df['Processes'].min(), df['Processes'].max(), 100)
        y_values = 1 / x_values
        plt.plot(x_values, y_values, linestyle='--', label='Linear Ref', color='black')

        plt.title('Scaling MPI and OMP versions')
        plt.xlabel('Proc/Threads')
        plt.ylabel('Times')
        plt.legend()
        plt.grid(True)

        # Salvare il grafico come PNG
        output_file = os.path.join(output_dir, f'confronto_scaling_gruppo_{i//2 + 1}.png')
        plt.savefig(output_file, bbox_inches='tight')
        print(f"Grafico salvato in: {output_file}")

        plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Utilizzo: python grafico_scaling_dual.py <path_al_file_csv>")
        sys.exit(1)

    file_path = sys.argv[1]
    plot_scaling_trends(file_path)
