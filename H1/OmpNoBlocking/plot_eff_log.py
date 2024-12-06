import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def generate_plots(input_file):

    data = pd.read_csv(input_file)
    

    columns_to_plot = ['SpeedupSym', 'SpeedupMatTrans', 'EfficiencySym', 'EfficiencyMatTrans']


    for dim in data['Dim'].unique():

        dim_data = data[data['Dim'] == dim]
        

        plt.figure(figsize=(8, 6))


        for column in columns_to_plot:
            plt.plot(dim_data['Threads'], dim_data[column], label=column, marker='o', linestyle='-', linewidth=2)
        

        max_threads = dim_data['Threads'].max()
        diagonal_x = range(0, max_threads + 1, 20)  
        diagonal_y = diagonal_x 
        plt.plot(diagonal_x, diagonal_y, label='Diagonal (0,0) (20,20) (40,40)...', linestyle='--', color='gray', linewidth=1)
        

        plt.yscale('log')
        

        plt.title(f"Comparison of {', '.join(columns_to_plot)} vs Threads for Dim {dim} (Log Scale)")
        plt.xlabel('Threads')
        plt.ylabel('Log of Value')
        plt.legend()
        plt.grid(True)


        output_plot_file = f"{input_file.split('.')[0]}_Dim{dim}_all_columns_with_diagonal_plot_log.png"
        plt.savefig(output_plot_file)
        plt.close()
        
        print(f"Combined plot for Dim {dim} saved as {output_plot_file}")

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python script.py <input_csv_file>")
        sys.exit(1)


    input_file = sys.argv[1]


    if not os.path.isfile(input_file):
        print(f"Error: File {input_file} does not exist.")
        sys.exit(1)


    generate_plots(input_file)

