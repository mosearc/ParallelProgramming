# Creazione dei grafici TRANS
for dim in unique_dims:
    subset = data[data['Dim'] == dim]
    plt.figure(figsize=(8, 6))
    
    # Traccia per la versione OMP (colore distintivo)
    omp_data = subset[subset['Version'] == 'OMP']
    if not omp_data.empty:
        plt.plot(omp_data['Threads'], omp_data['Trans'], marker='o', color='blue', label='OMP')
    
    # Linea costante per la versione SEQ (colore distintivo)
    seq_data = subset[subset['Version'] == 'SEQ']
    if not seq_data.empty and not seq_data['Trans'].isna().all():
        mean_seq_trans = seq_data['Trans'].mean()
        plt.axhline(y=mean_seq_trans, linestyle='--', color='green', label='SEQ (mean)')
    
    # Linea costante per la versione IMP (colore distintivo)
    imp_data = subset[subset['Version'] == 'IMP']
    if not imp_data.empty and not imp_data['Trans'].isna().all():
        mean_imp_trans = imp_data['Trans'].mean()
        plt.axhline(y=mean_imp_trans, linestyle='--', color='orange', label='IMP (mean)')
    
    # Impostazioni del grafico
    plt.title(f'Trans vs Threads per Dim = {dim}')
    plt.xlabel('Threads')
    plt.ylabel('Trans')
    plt.legend()
    plt.grid(True)
    plt.show()
