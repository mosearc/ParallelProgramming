import pandas as pd
import matplotlib.pyplot as plt
import sys

#python3 Speedup_Plot.py outputMeansLast_SpeedEff.csv 

def create_speedup_plots(csv_file):
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Get unique dimensions
    dimensions = df['Dim'].unique()
    
    # Create plots for each dimension
    for dim in dimensions:
        # Filter data for current dimension
        dim_data = df[df['Dim'] == dim]
        
        # Create plots for both SpeedupSym and SpeedupTrans
        for speedup_type in ['Sym', 'Trans']:
            plt.figure(figsize=(10, 6))
            
            # Get the speedup column name
            speedup_col = f'Speedup{speedup_type}'
            
            plt.plot(dim_data['Threads'], dim_data[speedup_col], 'b-', marker='o')
            
            # Customize the plot
            plt.title(f'Speedup {speedup_type} vs Threads (Dim = {dim})')
            plt.xlabel('Threads')
            plt.ylabel(f'Speedup {speedup_type}')
            plt.grid(True)
            
            # Add reference line y=x for ideal speedup
            max_threads = dim_data['Threads'].max()
            max_speedup = dim_data[speedup_col].max()
            ideal_line = min(max_threads, max_speedup)
            plt.plot([0, ideal_line], [0, ideal_line], 'r--', label='Ideal Speedup')
            
            plt.legend()
            
            # Save the plot
            plt.savefig(f'speedup_{speedup_type.lower()}_dim_{dim}.png')
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
