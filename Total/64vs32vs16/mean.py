import pandas as pd
import sys
import os

def process_csv(input_file):
    # Load the input file
    data = pd.read_csv(input_file)
    
    # Create new columns to store rolling means
    data['RollingMean_CheckSym'] = None
    data['RollingMean_MatTranspose'] = None

    # Calculate the mean every 5 rows for "CheckSym"
    for i in range(0, len(data), 5):
        mean_value_checksym = data['CheckSym'].iloc[i:i+5].mean()
        if pd.notna(mean_value_checksym):
            data.loc[i, 'RollingMean_CheckSym'] = mean_value_checksym
    
    # Calculate the mean every 5 rows for "MatTranspose"
    for i in range(0, len(data), 5):
        mean_value_mattranspose = data['MatTranspose'].iloc[i:i+5].mean()
        if pd.notna(mean_value_mattranspose):
            data.loc[i, 'RollingMean_MatTranspose'] = mean_value_mattranspose
    
    # Remove the original columns
    data = data.drop(columns=['CheckSym', 'MatTranspose'])
    
    # Remove rows where 'RollingMean_CheckSym' is None or NaN
    data = data[data['RollingMean_CheckSym'].notna()]
    
    # Create output file name
    base_name, ext = os.path.splitext(input_file)
    output_file = f"{base_name}_processed{ext}"
    
    # Save the cleaned data to the output file
    data.to_csv(output_file, index=False)
    print(f"Processed file saved as: {output_file}")

if __name__ == "__main__":
    # Check if the script is called with the right arguments
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)
    
    # Get input file path from arguments
    input_file = sys.argv[1]
    
    # Process the file
    process_csv(input_file)

