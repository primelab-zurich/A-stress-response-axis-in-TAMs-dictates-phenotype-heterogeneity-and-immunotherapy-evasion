import numpy as np
def gsea_dotplot(decoupler_res, fdr_cutoff = 0.2, title = 'GSEA', filename = None):
    import matplotlib.pyplot as plt
    import seaborn as sns
    fdr_label = 'pval'
    nes_label = 'norm'
    decoupler_res = decoupler_res[decoupler_res[fdr_label] <= fdr_cutoff]
    decoupler_res = decoupler_res[decoupler_res[nes_label] >= 0]
    sorted_dc_res = decoupler_res.sort_values(nes_label, ascending = True).reset_index(level=0)
    sorted_dc_res = sorted_dc_res[['source',nes_label, fdr_label]]
    
    if len(sorted_dc_res) <= 0:
        print('no Terms with positive normalized enrichment score')
        return
    max_nes = np.nanmax(sorted_dc_res[nes_label][sorted_dc_res[nes_label] != np.inf])
    min_nes = np.nanmin(sorted_dc_res[nes_label][sorted_dc_res[nes_label] != -np.inf])
    sorted_dc_res[nes_label][sorted_dc_res[nes_label] == np.inf] = max_nes
    sorted_dc_res[nes_label][sorted_dc_res[nes_label]== -np.inf] = min_nes

    color_scale_cutoff = fdr_cutoff*1.25
    norm = plt.Normalize(0, color_scale_cutoff)
    sm = plt.cm.ScalarMappable(cmap="Blues_r", norm=norm)
    sm.set_array([])
    
    fig, (ax) = plt.subplots(1, 1, figsize=(max_nes + 0.5, (len(sorted_dc_res))*0.24 + 1.9), dpi=150)
    fig.suptitle(title)
    ax = plt.scatter(
                    x=sorted_dc_res[nes_label],
                    y=sorted_dc_res['source'].astype(str) + ' '+ round(sorted_dc_res[fdr_label]*100).astype(np.int).astype(np.str_),
                    c=sorted_dc_res[fdr_label],
                    cmap='Blues_r',
                    norm = norm,
                    s=round(sorted_dc_res[nes_label]*240/(max_nes))
    )
    delta_nes = max_nes - min_nes
    plt.xlim([0, max_nes + 0.5])
    plt.ylim([-0.5, len(sorted_dc_res)-0.5])
    ax.figure.colorbar(sm)
    plt.tight_layout()
    if filename is not None:
        ax.figure.savefig(filename, transparent=False, bbox_inches = 'tight')

def gsea_decoupler_dotplot(gsea_results, pval_cutoff, title = 'GSEA', filename = None):
    import seaborn
    import numpy
    
    # fix index s.t. unused categories are removed and order is set acc. to norm
    gsea_results.reset_index(level = 'source', inplace=True)
    gsea_results.source.cat.remove_unused_categories(inplace=True)
    gsea_results = gsea_results.sort_values('norm', ascending = False)
    gsea_results = gsea_results.set_index('source')
    gsea_results.index = gsea_results.index.as_ordered()
    #print(gsea_results.index)
    gsea_results_filt = gsea_results[gsea_results.pval <= pval_cutoff]
    gsea_results_filt = gsea_results_filt.assign(
                **{"-log10(pval)": lambda x: -numpy.log10(x["pval"])}
            )
    p = seaborn.objects.Plot(
        data=(
            gsea_results_filt
        ),
        x="norm",
        y="source",
        color = "pval"
    ).add(seaborn.objects.Dot(pointsize=20)).scale(y = seaborn.objects.Nominal(order = gsea_results_filt.index), color = "Blues_r").save(filename)
    return p

from pathlib import Path
import pandas
def gmt_to_decoupler(pth: Path) -> pandas.DataFrame:
    """
    Parse a gmt file to a decoupler pathway dataframe.
    """
    from itertools import chain, repeat

    pathways = {}

    with Path(pth).open("r") as f:
        for line in f:
            name, _, *genes = line.strip().split("\t")
            pathways[name] = genes

    db =  pandas.DataFrame.from_records(
        chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
        columns=["geneset", "genesymbol"],
    )
    db.geneset = db['geneset'].astype('category')
    return db

def prepare_db_for_gsea(database, min_geneset_size = 15, max_geneset_size = 500):
    # Filtering genesets to match behaviour of fgsea
    geneset_size = database.groupby("geneset").size()
    gsea_genesets = geneset_size.index[(geneset_size > min_geneset_size) & (geneset_size < max_geneset_size)]
    database_filtered = database[~database.duplicated(("geneset", "genesymbol"))]
    database_filtered = database_filtered[database_filtered["geneset"].isin(gsea_genesets)]
    return database_filtered

import decoupler
def get_gsea_result(ranked_genes, database):
    score, norm, pval = decoupler.run_gsea(
    pandas.DataFrame(ranked_genes).T,
    database,
    source="geneset",
    target="genesymbol",
    min_n = 0
    )

    gsea_result = (
        pandas.concat({"score": score.T, "norm": norm.T, "pval": pval.T}, axis=1)
        .droplevel(level=1, axis=1)
    )
    return gsea_result