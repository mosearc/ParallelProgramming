import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def generate_plots(input_file):
    # Load the dataset
    data = pd.read_csv(input_file)
    
    # For each unique 'Dim', create a plot of RollingMean_CheckSym and RollingMean_MatTranspose vs Threads
    for dim in data['Dim'].unique():
        # Filter data for the current 'Dim'
        dim_data = data[data['Dim'] == dim]

        # Create a plot
        plt.figure(figsize=(8, 6))

        # Plot RollingMean_CheckSym vs Threads
        plt.plot(dim_data['Threads'], dim_data['RollingMean_CheckSym'], label='RollingMean_CheckSym', marker='o', linestyle='-', color='b')
        
        # Plot RollingMean_MatTranspose vs Threads
        plt.plot(dim_data['Threads'], dim_data['RollingMean_MatTranspose'], label='RollingMean_MatTranspose', marker='x', linestyle='--', color='r')
        
        # Add titles and labels
        plt.title(f"Rolling Means vs Threads for Dim {dim}")
        plt.xlabel('Threads')
        plt.ylabel('Value')
        plt.legend()
        plt.grid(True)

        # Save the plot as a PNG file
        output_plot_file = f"{input_file.split('.')[0]}_Dim{dim}_plot.png"
        plt.savefig(output_plot_file)
        plt.close()
        
        print(f"Plot for Dim {dim} saved as {output_plot_file}")

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

