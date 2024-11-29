import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def generate_plots(input_file):
    # Load the dataset
    data = pd.read_csv(input_file)
    
    # Filter SEQ version for reference lines
    seq_data = data[data['Version'] == 'SEQ']

    # For each unique 'Dim', create a plot of RollingMean_CheckSym and RollingMean_MatTranspose vs Threads
    for dim in data['Dim'].unique():
        # Filter data for the current 'Dim'
        dim_data = data[data['Dim'] == dim]

        # Get reference values for SEQ
        seq_ref = seq_data[seq_data['Dim'] == dim]
        seq_check_sym = seq_ref['RollingMean_CheckSym'].iloc[0] if not seq_ref.empty else None
        seq_mat_transpose = seq_ref['RollingMean_MatTranspose'].iloc[0] if not seq_ref.empty else None

        # Create a plot
        plt.figure(figsize=(8, 6))

        # Plot RollingMean_CheckSym vs Threads
        plt.plot(dim_data['Threads'], dim_data['RollingMean_CheckSym'], label='RollingMean_CheckSym', marker='o', linestyle='-', color='b')
        
        # Plot RollingMean_MatTranspose vs Threads
        plt.plot(dim_data['Threads'], dim_data['RollingMean_MatTranspose'], label='RollingMean_MatTranspose', marker='x', linestyle='--', color='r')
        
        # Add horizontal lines for SEQ reference values
        if seq_check_sym is not None:
            plt.axhline(y=seq_check_sym, color='green', linestyle='-', label='SEQ RollingMean_CheckSym')
        if seq_mat_transpose is not None:
            plt.axhline(y=seq_mat_transpose, color='orange', linestyle='-', label='SEQ RollingMean_MatTranspose')
        
        # Set logarithmic scale for Y-axis
        plt.yscale('log')

        # Add titles and labels
        plt.title(f"Rolling Means vs Threads for Dim {dim} (Log Scale)")
        plt.xlabel('Threads')
        plt.ylabel('Value (Log Scale)')
        plt.legend()
        plt.grid(True, which="both", linestyle='--', linewidth=0.5)

        # Save the plot as a PNG file
        output_plot_file = f"{input_file.split('.')[0]}_Dim{dim}_log_plot.png"
        plt.savefig(output_plot_file)
        plt.close()
        
        print(f"Log-scale plot for Dim {dim} saved as {output_plot_file}")

if __name__ == "__main__":
    # Check if the script is called with the correct arguments
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_csv_file>")
        sys.exit(1)

    # Get the input file from command line arguments
    input_file = sys.argv[1]

    # Validate the input file
    if not os.path.isfile(input_file):
        print(f"Error: File {input_file} does not exist.")
        sys.exit(1)

    # Generate plots
    generate_plots(input_file)

