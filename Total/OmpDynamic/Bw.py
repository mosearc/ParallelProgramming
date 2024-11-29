import matplotlib.pyplot as plt
import pandas as pd

# Carica i dati dal file CSV
# Sostituisci 'data.csv' con il nome del tuo file CSV
df = pd.read_csv('BW.csv')

# Estrai le colonne dal DataFrame
dim = df['Dim'].values
omp = df['Omp'].values
optimized = df['Optimized'].values
seq = df['Seq'].values
peak_value = df['Peak'].values[0]  # Assumiamo che Peak sia un valore costante

# Impostare la posizione per le barre (distribuite uniformemente)
x_positions = range(len(dim))  # Indici per l'asse X

# Creare il grafico a barre
plt.figure(figsize=(10, 6))
plt.bar([x - 0.2 for x in x_positions], omp, width=0.2, label='Omp', align='center')
plt.bar(x_positions, optimized, width=0.2, label='Optimized', align='center')
plt.bar([x + 0.2 for x in x_positions], seq, width=0.2, label='Seq', align='center')

# Aggiungere la linea orizzontale per il valore di Peak
if peak_value:
    plt.axhline(y=peak_value, color='red', linestyle='--', label=f'Peak = {peak_value:.2e}')

# Personalizzare il grafico
plt.xlabel('Dim')
plt.ylabel('Values')
plt.title('Comparison of Omp, Optimized, and Seq with Uniform Dim Distribution')
plt.xticks(x_positions, dim)
plt.legend()
plt.tight_layout()

# Salvare il grafico come immagine
plt.savefig('bar_plot_comparison.png')

# Mostrare il grafico
plt.show()

