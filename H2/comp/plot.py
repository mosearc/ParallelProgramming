import pandas as pd
import matplotlib.pyplot as plt
import sys


def generate_graphs(file_path):

    data = pd.read_csv(file_path)

    unique_dims = data['Dim'].unique()

    seq_data = data[data['Version'] == 'SEQ']

    for dim in unique_dims:
        plt.figure(figsize=(10, 6))

        seq_dim_data = seq_data[seq_data['Dim'] == dim]

        if not seq_dim_data.empty:
            plt.axhline(
                y=seq_dim_data['RollingMean_CheckSym'].iloc[0],
                color='blue', linestyle='--',
                label='SEQ - CheckSym'
            )
            plt.axhline(
                y=seq_dim_data['RollingMean_MatTranspose'].iloc[0],
                color='green', linestyle='--',
                label='SEQ - MatTranspose'
            )

        non_seq_data = data[(data['Dim'] == dim) & (data['Version'] != 'SEQ')]
        for version in non_seq_data['Version'].unique():
            version_data = non_seq_data[non_seq_data['Version'] == version]

            plt.plot(
                version_data['Threads/Proc'],
                version_data['RollingMean_CheckSym'],
                marker='o', linestyle='-', label=f'{version} - CheckSym'
            )
            plt.plot(
                version_data['Threads/Proc'],
                version_data['RollingMean_MatTranspose'],
                marker='s', linestyle='-', label=f'{version} - MatTranspose'
            )

        plt.title(f"Time for Dim={dim}")
        plt.xlabel("Threads/Proc")
        plt.ylabel("Time")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        #plt.show()

        output_path = f"Dim_{dim}_trend.png"
        plt.savefig(output_path)
        plt.close()
        print(f"Saved graph for Dim={dim} to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_csv_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    generate_graphs(file_path)
