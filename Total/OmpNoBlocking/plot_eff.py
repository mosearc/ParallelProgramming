import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def generate_plots(input_file):
    # Load the dataset
    data = pd.read_csv(input_file)
    
    # List of columns for which to generate plots
    columns_to_plot = ['SpeedupSym', 'SpeedupMatTrans', 'EfficiencySym', 'EfficiencyMatTrans']

    # For each unique 'Dim', create a combined plot for all the specified columns
    for dim in data['Dim'].unique():
        # Filter data for the current 'Dim'
        dim_data = data[data['Dim'] == dim]
        
        # Create a plot for all columns at once
        plt.figure(figsize=(8, 6))

        # Plot each of the specified columns vs Threads
        for column in columns_to_plot:
            plt.plot(dim_data['Threads'], dim_data[column], label=column, marker='o', linestyle='-', linewidth=2)
        
        # Plot the background diagonal line (0,0), (20,20), (40,40)...
        max_threads = dim_data['Threads'].max()
        diagonal_x = range(0, max_threads + 1, 20)  # Generate values from 0 to max_threads in steps of 20
        diagonal_y = diagonal_x  # y = x for the diagonal line
        plt.plot(diagonal_x, diagonal_y, label='Diagonal (0,0) (20,20) (40,40)...', linestyle='--', color='gray', linewidth=1)
        
        # Set the y-axis to a logarithmic scale
        #plt.yscale('log')
        
        # Add titles and labels
        plt.title(f"Comparison of {', '.join(columns_to_plot)} vs Threads for Dim {dim} (Log Scale)")
        plt.xlabel('Threads')
        plt.ylabel('Log of Value')
        plt.legend()
        plt.grid(True)

        # Save the plot as a PNG file
        output_plot_file = f"{input_file.split('.')[0]}_Dim{dim}_all_columns_with_diagonal_plot_log.png"
        plt.savefig(output_plot_file)
        plt.close()
        
        print(f"Combined plot for Dim {dim} saved as {output_plot_file}")

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

