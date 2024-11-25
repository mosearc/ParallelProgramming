import pandas as pd
import matplotlib.pyplot as plt
import sys

def create_speedup_plots(csv_file):
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Get unique dimensions
    dimensions = df['Dim'].unique()
    
    # Create a plot for each dimension
    for dim in dimensions:
        # Filter data for current dimension
        dim_data = df[df['Dim'] == dim]
        
        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.plot(dim_data['Threads'], dim_data['SpeedupSym'], 'b-', marker='o')
        
        # Customize the plot
        plt.title(f'Speedup Sym vs Threads (Dim = {dim})')
        plt.xlabel('Threads')
        plt.ylabel('Speedup Sym')
        plt.grid(True)
        
        # Add reference line y=x for ideal speedup
        max_threads = dim_data['Threads'].max()
        max_speedup = dim_data['SpeedupSym'].max()
        ideal_line = min(max_threads, max_speedup)
        plt.plot([0, ideal_line], [0, ideal_line], 'r--', label='Ideal Speedup')
        
        plt.legend()
        
        # Save the plot
        plt.savefig(f'speedup_dim_{dim}.png')
        plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <csv_file>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    try:
        create_speedup_plots(csv_file)
        print(f"Plots have been saved successfully!")
    except Exception as e:
        print(f"An error occurred: {e}")
