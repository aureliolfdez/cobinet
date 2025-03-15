import numpy as np
from sklearn.utils.validation import check_array #pip intall scikit-learn
from itertools import combinations #pip install itertools
import numpy as np #pip install numpy
import pandas as pd #pip install pandas
import os

def bcca(data, correlation_threshold=0.9, min_cols=3, dataset="", debug=False):

    genes = data.iloc[:, 0].astype(str).to_numpy()
    condiciones = np.array(list(map(str, data.columns[1:])))
    data = data.iloc[:, 1:].to_numpy()

    data = check_array(data, dtype=np.double, copy=True)
    
    num_rows, num_cols = data.shape
    biclusters = []

    for i, j in combinations(range(num_rows), 2):
        cols, corr = find_cols(data[i], data[j], correlation_threshold, min_cols)

        if len(cols) >= min_cols and corr >= correlation_threshold:
            rows = [i, j]

            for k, r in enumerate(data):
                if k != i and k != j and accept(data, rows, cols, r, correlation_threshold):
                    rows.append(k)

            bTemp = {
                    "genes": ",".join(sorted([genes[idx] for idx in rows])),
                    "condiciones": ",".join(sorted([condiciones[idx] for idx in cols]))
                    }

            if not exists(biclusters, bTemp):
                biclusters.append({"genes": ",".join([str(genes[idx]) for idx in rows]), "condiciones": ",".join([str(condiciones[idx]) for idx in cols])})
            
    saveResults(biclusters, genes, condiciones, data, dataset, correlation_threshold, min_cols, debug)   
    
    return biclusters

def find_cols(ri, rj, correlation_threshold, min_cols):
    
    cols = np.arange(len(ri), dtype=int)
    corr = corre(ri, rj)

    while corr < correlation_threshold and len(cols) >= min_cols:
        imax = find_max_decrease(ri, rj, cols)
        cols = np.delete(cols, imax)
        corr = corre(ri[cols], rj[cols])
    

    return cols, corr

def find_max_decrease(ri, rj, indices):
    
    kmax, greater = -1, float('-inf')

    for k in range(len(indices)):
        ind = np.concatenate((indices[:k], indices[k+1:]))
        result = corre(ri[ind], rj[ind])

        if result > greater:
            kmax, greater = k, result

    return kmax

def accept(data, rows, cols, r, correlation_threshold):
    
    for i in rows:
        corr = corre(r[cols], data[i, cols])

        if corr < correlation_threshold:
            return False

    return True

def corre(v, w):
    
    vc = v - np.mean(v)
    wc = w - np.mean(w)
    
    x = np.sum(vc * wc)
    y = np.sum(vc * vc) * np.sum(wc * wc)
    
    return np.abs(x / np.sqrt(y))

def exists(biclusters, bTemp):
    genes_sorted = ",".join(sorted(set(bTemp["genes"].split(","))))
    condiciones_sorted = ",".join(sorted(set(bTemp["condiciones"].split(","))))

    for b in biclusters:
        genes_b = ",".join(sorted(set(b["genes"].split(","))))
        condiciones_b = ",".join(sorted(set(b["condiciones"].split(","))))

        if genes_b == genes_sorted and condiciones_b == condiciones_sorted:
            return True  # Ya existe el bicluster en la lista

    return False  # No existe, se puede agregar

def saveResults(biclusters, genes, condiciones, data, dataset, correlation_threshold, min_cols, debug=False):

    output_dir = "results/bcca/" 
    output_file = os.path.join(output_dir, f"{dataset}_biclusters_cor-{correlation_threshold}_cols-{min_cols}.csv")

    biclusters_df = pd.DataFrame(biclusters)
    biclusters_df.to_csv(output_file, index=False, sep=";")
    print("Resultados guardados en 'results_biclusters.csv'", flush=True)

    biclusters_df = pd.read_csv(output_file, sep=';')
    grouped_biclusters = biclusters_df.groupby("condiciones")["genes"].apply(lambda x: ','.join(x)).reset_index()

    with open(f"results/bcca/{dataset}_infoNetworks_cor-{correlation_threshold}_cols-{min_cols}.txt", "w") as f:
        for index, row in grouped_biclusters.iterrows():
            f.write(f"{dataset}_network_{index}_cor-{correlation_threshold}_cols-{min_cols}.csv -> Columnas asociadas: {row['condiciones']}\n")
    print(f"Informacion de redes guardada en {dataset}_infoNetworks_cor-{correlation_threshold}_cols-{min_cols}.txt", flush=True)
            
    data_df = pd.DataFrame(data, index=genes, columns=condiciones)
    id=0
    for index, row in grouped_biclusters.iterrows():
        genes_group = list(set(row["genes"].split(',')))
        condiciones_group = row["condiciones"].split(',')
        subset_data = data_df.loc[genes_group, condiciones_group]
        correlations = []
        
        for g1, g2 in combinations(genes_group, 2):
            corr_value = np.corrcoef(subset_data.loc[g1], subset_data.loc[g2])[0, 1]
            correlations.append({"Interaction": g1, "name": g2, "Pearson": corr_value, "selected": "false", "shared interaction": g1, "shared name": g2, "Weight": corr_value})
        
        corr_df = pd.DataFrame(correlations)
        filename = f"results/bcca/{dataset}_network_{id}_cor-{correlation_threshold}_cols-{min_cols}.csv"
        corr_df.to_csv(filename, index=False, sep=',')
        print(f"Correlaciones guardadas en '{filename}'", flush=True)
        id+=1
     
    
    network_files = [f for f in os.listdir('./results/bcca') if f.startswith(f'{dataset}_network') and f.endswith('.csv')]
    
    if debug:
        print(network_files, flush=True)

    for network_file in network_files:
        
        if debug:
            print(network_file, flush=True) 
                 
        file_path = os.path.join('./results/bcca', network_file) 
        network_df = pd.read_csv(file_path)
        unique_genes = set(network_df['Interaction'].astype(str).tolist() + network_df['name'].astype(str).tolist())  
        genes_df = pd.DataFrame(unique_genes, columns=['Gene'])
        genes_file = file_path.replace('network', 'genes')
        genes_df.to_csv(genes_file, index=False, header=False)
        print(f'Archivo {genes_file} generado con éxito.', flush=True)