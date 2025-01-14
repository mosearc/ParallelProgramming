import pandas as pd
import argparse
import os

def calculate_scaling(input_file, output_file):
    df = pd.read_csv(input_file)


    required_columns = ['Processes', 'RollingMean_MatTranspose']
    if not all(col in df.columns for col in required_columns):
        raise ValueError("Il file deve contenere le colonne 'Processes' e 'RollingMean_MatTranspose'.")


    scaling = []
    scaling_value = None

    for index, row in df.iterrows():
        if row['Processes'] == 1:

            scaling_value = row['RollingMean_MatTranspose']
        if scaling_value is not None:
            scaling.append(scaling_value / row['RollingMean_MatTranspose'])
        else:
            scaling.append(None)

    df['scaling_RollingMean_MatTranspose'] = scaling

    df.to_csv(output_file, index=False)
    print(f"File salvato con successo come: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calcola la colonna "scaling" in un file CSV.')
    parser.add_argument('input_file', type=str, help='Percorso del file CSV di input')
    parser.add_argument('--output', type=str, default=None, help='Percorso del file CSV di output (opzionale)')

    args = parser.parse_args()


    output_file = args.output if args.output else os.path.splitext(args.input_file)[0] + '_with_scaling.csv'

    calculate_scaling(args.input_file, output_file)
