"""Annotate the 13 myeloid subtypes (microglia and perivascular macrophages) with 100 weighted markers per subtype.

    import annotate_myeloid_subtypes
    annotate_myeloid_subtypes(adata, weights='weighted_markers_top100.csv', batch='Source')

adds to adata, in place:

    adata.obs['pred_subtype']    one of the 13 subtypes
    adata.obs['pred_subclass']   Homeo, Adapt, ADAM, Prolif, exAM or PVM
    adata.obs['pred_prob']       probability of the predicted subtype
    adata.obsm['pred_log_prob']  log-probability of each of the 13 subtypes (cells x subtypes)

The same weights are used for cells (scRNA-seq) and nuclei (snRNA-seq). Needs anndata, numpy, pandas and scipy.
"""
import sys
import types
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.special import log_softmax

__all__ = ['annotate_myeloid_subtypes']

_HERE = Path(__file__).resolve().parent
_GENE_ID_COLUMNS = ('gene_id', 'gene_ids', 'featureid', 'ensembl_id')


def _raw_counts(adata, layer):
    """Raw counts (cells x genes, CSR) and the matching var table."""
    if layer is not None:
        counts, var = adata.layers[layer], adata.var
    elif 'counts' in adata.layers:
        counts, var = adata.layers['counts'], adata.var
    elif adata.raw is not None:
        counts, var = adata.raw.X, adata.raw.var
    else:
        counts, var = adata.X, adata.var
    counts = sp.csr_matrix(counts)
    if not np.allclose(counts.data[:100_000] % 1, 0):
        raise ValueError('raw counts are required: put them in adata.layers["counts"], adata.raw or adata.X, or name the layer')
    return counts, var


def annotate_myeloid_subtypes(adata, weights='weighted_markers_top100.csv', batch=None, layer=None):
    """Predict the myeloid subtype of every cell or nucleus from weighted marker genes.

    Parameters
    ----------
    adata
        AnnData of myeloid cells or nuclei only, e.g. ``adata[adata.obs['subclass'] == 'Micro_PVM'].copy()``, from
        one technology, with raw counts of all genes (the library size is the total count of a cell).
        Genes are matched by Ensembl ID when ``adata.var`` has them (index, or a column named gene_id, gene_ids,
        featureid or ensembl_id), otherwise by symbol.
    weights
        Weight table: the first row holds the intercepts, every other row one gene, one column per subtype.
        A name that is not found in the working directory is looked up next to this file.
    batch
        Column of ``adata.obs`` naming the batch (brain bank, site, chemistry). Every gene is centred and scaled
        within each batch; without it, over all of ``adata``. Each batch should hold a mix of myeloid states.
    layer
        Layer with the raw counts. Default: ``layers['counts']``, then ``adata.raw``, then ``adata.X``.

    Returns
    -------
    None. Adds ``pred_subtype``, ``pred_subclass`` and ``pred_prob`` to ``adata.obs`` and ``pred_log_prob``
    (a DataFrame, one column per subtype) to ``adata.obsm``.
    """
    path = Path(weights)
    if not path.exists() and (_HERE / path).exists():
        path = _HERE / path
    table = pd.read_csv(path)
    subtypes = list(table.columns[2:])
    intercept, table = table.iloc[0, 2:].to_numpy(dtype=float), table.iloc[1:]

    counts, var = _raw_counts(adata, layer)
    id_column = next((c for c in _GENE_ID_COLUMNS if c in var), None)
    if id_column is not None or str(var.index[0]).startswith('ENSG'):
        names = pd.Index(var[id_column] if id_column is not None else var.index).astype(str).str.replace(r'\.\d+$', '', regex=True)
        key = 'gene_id'
    else:
        names, key = pd.Index(var.index).astype(str), 'gene'
    table = table.dropna(subset=[key]).drop_duplicates(key).set_index(key)
    columns = np.flatnonzero(names.isin(table.index) & ~names.duplicated())
    if len(columns) < 0.5 * len(table):
        raise ValueError(f'only {len(columns)} of the {len(table)} marker genes were found in adata (matched by {key})')
    if len(columns) < 0.9 * len(table):
        warnings.warn(f'{len(columns)} of the {len(table)} marker genes were found in adata (matched by {key})')

    # log1p(CPM) of the marker genes, centred and scaled within each batch
    library = np.asarray(counts.sum(axis=1)).ravel()
    X = np.log1p((sp.diags(1e6 / np.maximum(library, 1)) @ counts[:, columns]).toarray()).astype(np.float32)
    groups = adata.obs[batch].astype(str).to_numpy() if batch is not None else np.zeros(adata.n_obs, dtype=int)
    for group in np.unique(groups):
        rows = np.flatnonzero(groups == group)
        if len(rows) < 100:
            warnings.warn(f'batch {group!r} has only {len(rows)} cells: scaling within it is unreliable')
        sd = X[rows].std(axis=0)
        X[rows] = (X[rows] - X[rows].mean(axis=0)) / np.where(sd > 0, sd, 1)

    # weighted sum of the markers for each subtype -> probabilities
    W = table.loc[names[columns], subtypes].to_numpy(dtype=np.float32)
    log_prob = log_softmax(X @ W + intercept, axis=1).astype(np.float32)
    best = log_prob.argmax(axis=1)
    subclasses = list(dict.fromkeys(s.split('_')[0] for s in subtypes))
    adata.obs['pred_subclass'] = pd.Categorical([subtypes[i].split('_')[0] for i in best], categories=subclasses)
    adata.obs['pred_subtype'] = pd.Categorical.from_codes(best, categories=subtypes)
    adata.obs['pred_prob'] = np.exp(log_prob.max(axis=1))
    adata.obsm['pred_log_prob'] = pd.DataFrame(log_prob, index=adata.obs_names, columns=subtypes)


class _CallableModule(types.ModuleType):
    """Lets ``import annotate_myeloid_subtypes; annotate_myeloid_subtypes(adata)`` work."""

    def __call__(self, *args, **kwargs):
        return annotate_myeloid_subtypes(*args, **kwargs)


sys.modules[__name__].__class__ = _CallableModule
