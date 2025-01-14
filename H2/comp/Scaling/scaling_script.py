import pandas as pd
import argparse
import os

def calculate_scaling(input_file, output_file):

    df = pd.read_csv(input_file)

    required_columns = ['Processes', 'RollingMean_MatTranspose', 'RollingMean_CheckSym']
    if not all(col in df.columns for col in required_columns):
        raise ValueError("Il file deve contenere le colonne 'Processes', 'RollingMean_MatTranspose' e 'RollingMean_CheckSym'.")

    scaling_mat = []
    scaling_value_mat = None

    for index, row in df.iterrows():
        if row['Processes'] == 1:
            scaling_value_mat = row['RollingMean_MatTranspose']
        if scaling_value_mat is not None:
            scaling_mat.append(scaling_value_mat / row['RollingMean_MatTranspose'])
        else:
            scaling_mat.append(None)

    df['scaling_RollingMean_MatTranspose'] = scaling_mat


    scaling_sym = []
    scaling_value_sym = None

    for index, row in df.iterrows():
        if row['Processes'] == 1:
            scaling_value_sym = row['RollingMean_CheckSym']
        if scaling_value_sym is not None:
            scaling_sym.append(scaling_value_sym / row['RollingMean_CheckSym'])
        else:
            scaling_sym.append(None)

    df['scaling_RollingMean_CheckSym'] = scaling_sym

    df.to_csv(output_file, index=False)
    print(f"File salvato con successo come: {output_file}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calcola le colonne "scaling_RollingMean_MatTranspose" e "scaling_RollingMean_CheckSym" in un file CSV.')
    parser.add_argument('input_file', type=str, help='Percorso del file CSV di input')
    parser.add_argument('--output', type=str, default=None, help='Percorso del file CSV di output (opzionale)')

    args = parser.parse_args()

    output_file = args.output if args.output else os.path.splitext(args.input_file)[0] + '_with_scaling.csv'

    calculate_scaling(args.input_file, output_file)
