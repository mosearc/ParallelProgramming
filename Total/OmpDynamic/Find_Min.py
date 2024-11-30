import pandas as pd
import sys
import os

def find_min_rolling_mean_mattrans(input_file):

    data = pd.read_csv(input_file)


    omp_data = data[data['Version'] == 'OMP']


    min_values = omp_data.groupby('Dim')['RollingMean_MatTranspose'].min()


    print("Minimi valori di RollingMean_MatTranspose per ogni Dim (Version = OMP):")
    print(min_values)


    output_file = 'Min_MatTranspose.csv'
    min_values.to_csv(output_file)
    print(f"Risultati salvati in {output_file}")

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Uso: python script.py <file_csv_input>")
        sys.exit(1)


    input_file = sys.argv[1]


    if not os.path.isfile(input_file):
        print(f"Errore: Il file {input_file} non esiste.")
        sys.exit(1)


    find_min_rolling_mean_mattrans(input_file)

