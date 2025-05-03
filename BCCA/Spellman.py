from BCCA import bcca
import pandas as pd
import numpy as np

############
# Spellman #
############
# Spellman: 24 cols
normalColsDataset = 24

#percentageCols = [1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]
#correlationThresholds = [0.9, 0.85, 0.8, 0.75, 0.7]

percentageCols = [0.70]
correlationThresholds = [0.7]

# Diccionario con los datasets
datasets = {
    "normal": "../datasets/Spellman.csv"
}

minGenes = 2        # Numero de genes minimo por cada red.
interactMin = 1 # Numero de interacciones mínimas por cada red.

# Procesamiento de cada combinación
for threshold in correlationThresholds:
    print(f"Processing with correlation_threshold = {threshold}", flush=True)
    
    for dataset_type, file_path in datasets.items():    
        for percentage in percentageCols:
            cols = round(normalColsDataset * percentage)
            print(f"{dataset_type.capitalize()} - {int(percentage)}% cols: {cols}", flush=True)
            
            data = pd.read_csv(file_path, sep=",")
            bcca(data, correlation_threshold=threshold, min_cols=cols, dataset="Spellman", geneMin=minGenes, interactionsMin=interactMin)
