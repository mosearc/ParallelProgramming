import pandas as pd
import matplotlib.pyplot as plt
import sys

# Function to generate graphs
def generate_graphs(file_path):
    # Load the dataset
    data = pd.read_csv(file_path)

    # Get unique 'Dim' values
    unique_dims = data['Dim'].unique()

    # Filter SEQ version data
    seq_data = data[data['Version'] == 'SEQ']

    # Plot for each unique Dim value
    for dim in unique_dims:
        plt.figure(figsize=(10, 6))

        # SEQ version for the current Dim
        seq_dim_data = seq_data[seq_data['Dim'] == dim]

        # Plot SEQ lines for RollingMean_CheckSym and RollingMean_MatTranspose
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

        # Non-SEQ versions for the same Dim
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

        # Set y-axis to logarithmic scale
        plt.yscale('log')

        # Plot formatting
        plt.title(f"Time for Dim={dim} (Log Scale)")
        plt.xlabel("Threads/Proc")
        plt.ylabel("Time (Log Scale)")
        plt.legend()
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.tight_layout()
        #plt.show()

        # Save the graph as a PNG file
        output_path = f"Dim_{dim}_trend.png"
        plt.savefig(output_path)
        plt.close()
        print(f"Saved graph for Dim={dim} to {output_path}")

# Main Execution: Get input file from command line
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_csv_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    generate_graphs(file_path)
