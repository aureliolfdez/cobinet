from BCCA import bcca
import pandas as pd
import numpy as np

########
# FC 1 #
########
# Normal DEGS: 13 cols
# Tumor DEGs: 44 cols
normalColsDataset = 18
tumorColsDataset = 18

percentageCols = [1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]
correlationThresholds = [0.9, 0.85, 0.8, 0.75, 0.7]

# Diccionario con los datasets
datasets = {
    "normal": "datasets/GSE99671/degs/NormalSarcoma_DEGs_table.csv",
    "tumor": "datasets/GSE99671/degs/tumorSarcoma_DEGs_table.csv"
}

minGenes = 10       # Numero de genes minimo por cada red.
interactMin = 10    # Numero de interacciones mínimas por cada red.

# Procesamiento de cada combinación
for threshold in correlationThresholds:
    print(f"Processing with correlation_threshold = {threshold}", flush=True)
    
    for dataset_type, file_path in datasets.items():
        cols_total = normalColsDataset if dataset_type == "normal" else tumorColsDataset
        
        for percentage in percentageCols:
            cols = round(cols_total * percentage)
            print(f"{dataset_type.capitalize()} - {int(percentage)}% cols: {cols}", flush=True)
            
            data = pd.read_csv(file_path, sep="\t")
            bcca(data, correlation_threshold=threshold, min_cols=cols, dataset="GSE99671", geneMin=minGenes, interactionsMin=interactMin)