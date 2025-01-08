import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def generate_plots(input_file):

    data = pd.read_csv(input_file)

    required_columns = ['Dim', 'Threads/Proc', 'Version', 'SpeedupSym', 'SpeedupMatTrans', 'EfficiencySym', 'EfficiencyMatTrans']
    if not all(col in data.columns for col in required_columns):
        raise ValueError("Il file deve contenere le colonne: " + ", ".join(required_columns))

    columns_to_plot = ['SpeedupSym', 'SpeedupMatTrans', 'EfficiencySym', 'EfficiencyMatTrans']

    for dim in data['Dim'].unique():
        dim_data = data[data['Dim'] == dim]

        omp_data = dim_data[dim_data['Version'] == 'OMP']
        mpi_data = dim_data[dim_data['Version'] == 'MPI DT']

        plt.figure(figsize=(12, 8))

        for column in columns_to_plot:
            if not omp_data.empty:
                plt.plot(omp_data['Threads/Proc'], omp_data[column], label=f'OMP - {column}', marker='o', linestyle='-', linewidth=2)
            if not mpi_data.empty:
                plt.plot(mpi_data['Threads/Proc'], mpi_data[column], label=f'MPI - {column}', marker='s', linestyle='--', linewidth=2)

        max_threads = dim_data['Threads/Proc'].max()
        diagonal_x = range(0, max_threads + 1, max(1, max_threads // 10))
        diagonal_y = diagonal_x
        plt.plot(diagonal_x, diagonal_y, label='Diagonal', linestyle='--', color='gray', linewidth=1)

        plt.yscale('log')

        plt.title(f"Speedup & Efficiency for Dim {dim}")
        plt.xlabel('Threads/Proc')
        plt.ylabel('Times')
        plt.legend()
        plt.grid(True)

        output_plot_file = f"{os.path.splitext(input_file)[0]}_Dim{dim}_OMP_MPI_comparison_plot.png"
        plt.savefig(output_plot_file)
        plt.close()

        print(f"Plot for Dim {dim} saved as {output_plot_file}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_dim_versions.py <input_csv_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    if not os.path.isfile(input_file):
        print(f"Error: File {input_file} does not exist.")
        sys.exit(1)

    generate_plots(input_file)
