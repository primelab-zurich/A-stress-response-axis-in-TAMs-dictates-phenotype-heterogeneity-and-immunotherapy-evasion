import scanpy as sc
import pathlib
import anndata
import pandas as pd
def read_visium_hd(sample_path: pathlib.Path) -> anndata.AnnData:
    tissue_positions_csv_path = sample_path / 'spatial' / 'tissue_positions_list.csv'
    if not tissue_positions_csv_path.exists():
        tissue_positions = pd.read_parquet(sample_path / 'spatial' / 'tissue_positions.parquet')
        tissue_positions.to_csv(tissue_positions_csv_path, index = False)

    adata = sc.read_visium(sample_path)
    adata.var_names_make_unique()
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(float)
    return adata

def preprocess(sample_path: pathlib.Path, preprocessed_path: pathlib.Path):
    ad = read_visium_hd(sample_path)
    sc.pp.normalize_total(ad)
    sc.pp.log1p(ad)
    ad.write(preprocessed_path)
    