import pandas as pd

#script per filtrare outputMean e renderlo outputMeanFinal  python3 filterOutput.py 

# Carica il file CSV
file_path = 'outputMeansLast.csv'  # Sostituisci con il percorso del tuo file
data = pd.read_csv(file_path)

# Filtra le righe con valori non nulli nella colonna 'Trans'
filtered_data = data.dropna(subset=['Trans'])

# Rimuovi le righe con "IMP" nella colonna 'Version'
filtered_data_no_imp = filtered_data[filtered_data['Version'] != 'IMP']

# Elimina le colonne 'CheckSym' e 'MatTranspose'
final_data = filtered_data_no_imp.drop(columns=['CheckSym', 'MatTranspose'])

# Salva il dataset aggiornato in un nuovo file CSV
output_path = 'outputMeansLast_final.csv'  # Sostituisci con il percorso desiderato per il file salvato
final_data.to_csv(output_path, index=False)

print(f"File filtrato salvato come {output_path}")

