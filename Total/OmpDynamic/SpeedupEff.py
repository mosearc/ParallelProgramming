import pandas as pd
import sys
import os

def process_csv(input_file, output_file):
    # Load the dataset
    data = pd.read_csv(input_file)

    # Add SpeedupSym column
    data['SpeedupSym'] = None
    for index, row in data.iterrows():
        if row['Version'] == 'OMP':
            seq_value = data.loc[(data['Dim'] == row['Dim']) & (data['Version'] == 'SEQ'), 'RollingMean_CheckSym']
            if not seq_value.empty and row['RollingMean_CheckSym'] != 0:
                data.at[index, 'SpeedupSym'] = seq_value.values[0] / row['RollingMean_CheckSym']

    # Add SpeedupMatTrans column
    data['SpeedupMatTrans'] = None
    for index, row in data.iterrows():
        if row['Version'] == 'OMP':
            seq_value = data.loc[(data['Dim'] == row['Dim']) & (data['Version'] == 'SEQ'), 'RollingMean_MatTranspose']
            if not seq_value.empty and row['RollingMean_MatTranspose'] != 0:
                data.at[index, 'SpeedupMatTrans'] = seq_value.values[0] / row['RollingMean_MatTranspose']

    # Add EfficiencySym column
    data['EfficiencySym'] = data.apply(
        lambda row: (float(row['SpeedupSym']) / row['Threads']) * 100 
        if row['SpeedupSym'] is not None else None, axis=1
    )

    # Add EfficiencyMatTrans column
    data['EfficiencyMatTrans'] = data.apply(
        lambda row: (float(row['SpeedupMatTrans']) / row['Threads']) * 100 
        if row['SpeedupMatTrans'] is not None else None, axis=1
    )

    # Save the result to a new CSV file
    data.to_csv(output_file, index=False)
    print(f"Processed data has been saved to {output_file}")

if __name__ == "__main__":
    # Check if the script is called with the correct arguments
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_csv_file>")
        sys.exit(1)

    # Get the input file from command line arguments
    input_file = sys.argv[1]

    # Validate the input file
    if not os.path.isfile(input_file):
        print(f"Error: File {input_file} does not exist.")
        sys.exit(1)

    # Create the output file name by adding "SpeeEff" before the extension
    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_SpeeEff{ext}"

    # Process the CSV
    process_csv(input_file, output_file)

