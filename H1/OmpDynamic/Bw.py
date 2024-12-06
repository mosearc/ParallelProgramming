import matplotlib.pyplot as plt
import pandas as pd


df = pd.read_csv('BW.csv')


dim = df['Dim'].values
omp = df['Omp'].values
optimized = df['Optimized'].values
seq = df['Seq'].values
peak_value = df['Peak'].values[0]  


x_positions = range(len(dim)) 


plt.figure(figsize=(10, 6))
plt.bar([x - 0.2 for x in x_positions], omp, width=0.2, label='Omp', align='center')
plt.bar(x_positions, optimized, width=0.2, label='Optimized', align='center')
plt.bar([x + 0.2 for x in x_positions], seq, width=0.2, label='Seq', align='center')


if peak_value:
    plt.axhline(y=peak_value, color='red', linestyle='--', label=f'Peak = {peak_value:.2e}')


plt.xlabel('Dim')
plt.ylabel('Values')
plt.title('Comparison of Omp, Optimized, and Seq with Uniform Dim Distribution')
plt.xticks(x_positions, dim)
plt.legend()
plt.tight_layout()


plt.savefig('bar_plot_comparison.png')


plt.show()

