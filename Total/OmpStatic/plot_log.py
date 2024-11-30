import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def generate_plots(input_file):

    data = pd.read_csv(input_file)
    

    seq_data = data[data['Version'] == 'SEQ']

 
    for dim in data['Dim'].unique():

        dim_data = data[data['Dim'] == dim]


        seq_ref = seq_data[seq_data['Dim'] == dim]
        seq_check_sym = seq_ref['RollingMean_CheckSym'].iloc[0] if not seq_ref.empty else None
        seq_mat_transpose = seq_ref['RollingMean_MatTranspose'].iloc[0] if not seq_ref.empty else None


        plt.figure(figsize=(8, 6))


        plt.plot(dim_data['Threads'], dim_data['RollingMean_CheckSym'], label='RollingMean_CheckSym', marker='o', linestyle='-', color='b')
        

        plt.plot(dim_data['Threads'], dim_data['RollingMean_MatTranspose'], label='RollingMean_MatTranspose', marker='x', linestyle='--', color='r')
        

        if seq_check_sym is not None:
            plt.axhline(y=seq_check_sym, color='green', linestyle='-', label='SEQ RollingMean_CheckSym')
        if seq_mat_transpose is not None:
            plt.axhline(y=seq_mat_transpose, color='orange', linestyle='-', label='SEQ RollingMean_MatTranspose')
        

        plt.yscale('log')


        plt.title(f"Rolling Means vs Threads for Dim {dim} (Log Scale)")
        plt.xlabel('Threads')
        plt.ylabel('Value (Log Scale)')
        plt.legend()
        plt.grid(True, which="both", linestyle='--', linewidth=0.5)


        output_plot_file = f"{input_file.split('.')[0]}_Dim{dim}_log_plot.png"
        plt.savefig(output_plot_file)
        plt.close()
        
        print(f"Log-scale plot for Dim {dim} saved as {output_plot_file}")

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python script.py <input_csv_file>")
        sys.exit(1)


    input_file = sys.argv[1]


    if not os.path.isfile(input_file):
        print(f"Error: File {input_file} does not exist.")
        sys.exit(1)


    generate_plots(input_file)

