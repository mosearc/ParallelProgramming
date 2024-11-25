import pandas as pd
import matplotlib.pyplot as plt
import sys

def create_efficiency_plots(csv_file):
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Get unique dimensions
    dimensions = df['Dim'].unique()
    
    # Create plots for each dimension
    for dim in dimensions:
        # Filter data for current dimension
        dim_data = df[df['Dim'] == dim]
        
        # Create plots for both EfficencySym and EfficencyTrans
        for eff_type in ['Sym', 'Trans']:
            plt.figure(figsize=(10, 6))
            
            # Get the efficiency column name
            eff_col = f'Efficency{eff_type}'
            
            plt.plot(dim_data['Threads'], dim_data[eff_col], 'b-', marker='o')
            
            # Customize the plot
            plt.title(f'Efficiency {eff_type} vs Threads (Dim = {dim})')
            plt.xlabel('Threads')
            plt.ylabel(f'Efficiency {eff_type}')
            plt.grid(True)
            
            # Set y-axis limits from 0 to max with some padding
            max_eff = dim_data[eff_col].max()
            plt.ylim(0, max_eff * 1.1)
            
            # Save the plot
            plt.savefig(f'efficiency_{eff_type.lower()}_dim_{dim}.png')
            plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <csv_file>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    try:
        create_efficiency_plots(csv_file)
        print(f"Efficiency plots have been saved successfully!")
    except Exception as e:
        print(f"An error occurred: {e}")
