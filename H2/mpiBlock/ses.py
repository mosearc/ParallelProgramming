import pandas as pd
import sys
import os

def process_csv(input_file, output_file):
    data = pd.read_csv(input_file)

    data['SpeedupMatTrans'] = None
    for index, row in data.iterrows():
        if row['Version'] in ['MPI BDV', 'MPI DT']:
            seq_value = data.loc[(data['Dim'] == row['Dim']) & (data['Version'] == 'SEQ'), 'RollingMean_MatTranspose']
            if not seq_value.empty and row['RollingMean_MatTranspose'] != 0:
                data.at[index, 'SpeedupMatTrans'] = seq_value.values[0] / row['RollingMean_MatTranspose']


    data['EfficiencyMatTrans'] = data.apply(
        lambda row: (float(row['SpeedupMatTrans']) / row['Processes']) * 100
        if row['SpeedupMatTrans'] is not None and row['Version'] in ['MPI BDV', 'MPI DT'] else None,
        axis=1
    )



    data.to_csv(output_file, index=False)
    print(f"Processed data has been saved to {output_file}")

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python script.py <input_csv_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    if not os.path.isfile(input_file):
        print(f"Error: File {input_file} does not exist.")
        sys.exit(1)

    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_SpeeEff{ext}"

    process_csv(input_file, output_file)
