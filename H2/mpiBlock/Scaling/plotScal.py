import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def plot_scaling_trends(file_path):

    df = pd.read_csv(file_path)

    output_dir = 'grafici_scaling'
    os.makedirs(output_dir, exist_ok=True)

    df['Version_Change'] = df['Version'] != df['Version'].shift()
    df['Group'] = df['Version_Change'].cumsum()

    version_groups = df.groupby('Group')
    grouped_versions = list(version_groups)

    for i in range(0, len(grouped_versions), 2):
        plt.figure(figsize=(12, 7))

        for j in range(2):
            if i + j < len(grouped_versions):
                group_id, group = grouped_versions[i + j]
                version = group['Version'].iloc[0]
                plt.plot(group['Processes'], group['scaling_RollingMean_MatTranspose'],
                         marker='o', label=f'{version} ')

        x_values = np.linspace(df['Processes'].min(), df['Processes'].max(), 100)
        y_values = 1 / x_values
        plt.plot(x_values, y_values, linestyle='--', label='Linear Ref', color='black')

        plt.title('Scaling MPI Block & MPI')
        plt.xlabel('Processes')
        plt.ylabel('Times')
        plt.legend()
        plt.grid(True)

        output_file = os.path.join(output_dir, f'confronto_scaling_gruppo_{i//2 + 1}.png')
        plt.savefig(output_file, bbox_inches='tight')
        print(f"Grafico salvato in: {output_file}")

        plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Utilizzo: python grafico_scaling.py <path_al_file_csv>")
        sys.exit(1)

    file_path = sys.argv[1]
    plot_scaling_trends(file_path)
