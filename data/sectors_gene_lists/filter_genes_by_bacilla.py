import pandas as pd
import os
from time import perf_counter

def read_sector_gene_list(sector_filename):
    """Reads a sector gene list file and returns a list of gene IDs."""
    t0 = perf_counter()
    df = pd.read_csv(sector_filename, header=None)
    elapsed = perf_counter() - t0
    print(f"Loaded {len(df)} genes from {sector_filename} (read in {elapsed:.3f}s)")
    print(df.head())
    gene_list = df[0].tolist()
    return gene_list

if __name__ == "__main__":
    total_start = perf_counter()

    # Load Bacillaria gene list
    t0 = perf_counter()
    bacilla_gene_list = pd.read_csv('../filter/ID_bacilla.csv', header=0).iloc[:,0].tolist()
    read_time = perf_counter() - t0
    print(f"Loaded {len(bacilla_gene_list)} Bacillaria genes in {read_time:.3f}s.")
    print(bacilla_gene_list[:5])

    # Prepare set once
    bacilla_set = set(bacilla_gene_list)

    # Define sector filenames
    sector_filenames = {
        'A': 'genes_in_sector_A.txt',
        'C': 'genes_in_sector_C.txt',
        'P': 'genes_in_sector_P.txt',
        'Q': 'genes_in_sector_Q.txt',
        'R': 'genes_in_sector_R.txt'
    }

    for sector_key, filename in sector_filenames.items():
        sector_start = perf_counter()
        sector_gene_list = read_sector_gene_list(filename)
        filtered_genes = [gene for gene in sector_gene_list if gene in bacilla_set]
        output_filename = f'filtered_{filename}'
        pd.DataFrame(filtered_genes).to_csv(output_filename, index=False, header=False)
        sector_elapsed = perf_counter() - sector_start
        print(f"Sector {sector_key}: filtered {len(filtered_genes)} genes -> {output_filename} (time {sector_elapsed:.3f}s)")

    total_elapsed = perf_counter() - total_start
    print(f"Total runtime: {total_elapsed:.3f}s")