import pandas as pd
import sys
import os

def process_csv(input_file):

    data = pd.read_csv(input_file)

    data['RollingMean_MatTranspose'] = None


    

    for i in range(0, len(data), 5):
        mean_value_mattranspose = data['MatTranspose'].iloc[i:i+5].mean()
        if pd.notna(mean_value_mattranspose):
            data.loc[i, 'RollingMean_MatTranspose'] = mean_value_mattranspose
    

    data = data.drop(columns=['CheckSym', 'MatTranspose'])
    

    data = data[data['RollingMean_MatTranspose'].notna()]
    

    base_name, ext = os.path.splitext(input_file)
    output_file = f"{base_name}_processed{ext}"

    data.to_csv(output_file, index=False)
    print(f"Processed file saved as: {output_file}")

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)
    

    input_file = sys.argv[1]
    

    process_csv(input_file)

