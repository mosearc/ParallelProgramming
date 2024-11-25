import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#python3 BmemPlotLog.py Bmem.csv

# Leggi il file CSV
df = pd.read_csv('Bmem.csv')

# Crea il grafico
plt.figure(figsize=(12, 8))

# Imposta la scala logaritmica per l'asse y
plt.yscale('log')

# Posizione delle barre
x = np.arange(len(df['dim']))
width = 0.25  # Larghezza delle barre

# Crea le barre per ogni metrica
plt.bar(x - width, df['omp'], width, label='OMP', color='skyblue')
plt.bar(x, df['imp'], width, label='IMP', color='lightgreen')
plt.bar(x + width, df['seq'], width, label='SEQ', color='salmon')

# Aggiungi la linea orizzontale per il peak value (prendi il primo valore non nullo)
peak_value = df['peak'].dropna().iloc[0]
plt.axhline(y=peak_value, color='red', linestyle='--', label=f'Peak ({peak_value:,.0f})')

# Personalizza il grafico
plt.xlabel('Dimensione')
plt.ylabel('Valore (scala log)')
plt.title('Confronto Performance')
plt.xticks(x, df['dim'], rotation=45)
plt.legend(loc='upper right')

# Formatta gli assi per numeri grandi
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:,.0f}'))

# Aggiusta il layout per evitare sovrapposizioni
plt.tight_layout()

# Salva il grafico come PNG
plt.savefig('performance_plot.png', dpi=300, bbox_inches='tight')
