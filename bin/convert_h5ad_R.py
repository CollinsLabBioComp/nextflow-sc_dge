#!/usr/bin/env python


__version__ = '0.0.1'


import argparse
import os
from packaging.version import Version
import scipy
import scipy.io
import gzip
import pandas as pd
import numpy as np
import scanpy as sc


def is_interaction(term):
    return(len(term.split(':')) > 1)


def convert_h5ad_R(
    adata,
    out_dir='tenx_metadata',
    verbose=True
):
    """Write 10x like data from h5ad data.

    Include a cell_metadata file.

    Parameters
    ----------
    adata : pandas.DataFrame
        Description of parameter `adata`.
    out_dir : string
        Description of parameter `out_dir`.
    verbose : boolean
        Description of parameter `verbose`.

    Returns
    -------
    execution_code : int
    """
    # Make the output directory if it does not exist.
    if out_dir == '':
        out_dir = os.getcwd()
    else:
        os.makedirs(out_dir, exist_ok=True)

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if Version(pd.__version__) > Version('1.0.0'):
        compression_opts = dict(method='gzip', compresslevel=9)

    # Save the barcodes.
    out_f = os.path.join(
        out_dir,
        'barcodes.tsv.gz'
    )
    if verbose:
        print('Writing {}'.format(out_f))
    pd.DataFrame(adata.obs.index).to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the features.
    # If column 'gene_symbols' exists, expect Ensembl IDs as index
    # If column 'gene_ids' exists, expect gene symbols as index
    out_f = os.path.join(
        out_dir,
        'features.tsv.gz'
    )
    if verbose:
        print('Writing {}'.format(out_f))

    # Check if read in by gene_symbol or gene_id
    if 'gene_symbols' in adata.var.columns:
        df_features = pd.DataFrame(data=[
            adata.var.index.values,
            adata.var.loc[:, 'gene_symbols'].values,
            ["Gene Expression"] * adata.n_vars
        ]).T
    elif 'gene_ids' in adata.var.columns:
        df_features = pd.DataFrame(data=[
            adata.var.loc[:, 'gene_ids'].values,
            adata.var.index.values,
            ["Gene Expression"] * adata.n_vars
        ]).T
    else:
        raise Exception(
            'Could not find "gene_symbols" or "gene_ids" in adata.var'
        )

    df_features.to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the  matrix
    out_mtx = adata.X.transpose()
    if not isinstance(out_mtx, scipy.sparse.csr_matrix):
        out_mtx = scipy.sparse.csr_matrix(out_mtx)
    out_f = os.path.join(
        out_dir,
        'matrix.mtx.gz'
    )
    if verbose:
        print('Writing {}'.format(out_f))
        # print(out_mtx.dtype)  # it looks like even counts stored as float
    with gzip.open(out_f, 'wb', compresslevel=9) as f:
        scipy.io.mmwrite(
            f,
            out_mtx,
            comment='metadata_json: {"format_version": 2}'
            # field='real'  # can be integer if counts otherwise real
        )

    if verbose:
        print('Done.')

    return 0


def write_metadata(
    adata,
    cols_to_retain,
    out_dir='tenx_metadata',
    verbose=True
):
    """Write h5ad metadata

    Parameters
    ----------
    adata : pandas.DataFrame
        Description of parameter `adata`.
    out_dir : string
        Description of parameter `out_dir`.
    verbose : boolean
        Description of parameter `verbose`.

    Returns
    -------
    execution_code : int
    """
    # Make the output directory if it does not exist.
    if out_dir == '':
        out_dir = os.getcwd()
    else:
        os.makedirs(out_dir, exist_ok=True)

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if Version(pd.__version__) > Version('1.0.0'):
        compression_opts = dict(method='gzip', compresslevel=9)

    # Now write cell metadata as TSV
    out_f = os.path.join(
        out_dir,
        'cell_metadata.tsv.gz'
    )
    if verbose:
        print('Writing {}'.format(out_f))

    cols_intersect = np.intersect1d(adata.obs.columns.tolist(), cols_to_retain)
    adata.obs[cols_intersect].to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=True,
        header=True
    )
    col_diff = np.setdiff1d(cols_intersect, cols_to_retain)
    if verbose and len(col_diff) > 0:
        print(
            ('Skipped the following columns because they were not '
             'in dataframe: {}').format(', '.join(col_diff))
        )

    if verbose:
        print('Done.')

    return 0


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read an h5ad object and write matrix files similar to 10x
            output. Additionally, write out cell and feature metadata
            to be read in R. This script requires pandas >1.0.1
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5ad', '--h5ad_file',
        action='store',
        dest='h5ad',
        required=True,
        help='H5AD file.'
    )

    parser.add_argument(
        '-target', '--variable_target',
        action='store',
        dest='target',
        required=True,
        help='Column to test.'
    )

    parser.add_argument(
        '-cl', '--cell_label',
        action='store',
        dest='cl',
        default='all',
        help=''' Comma-separated list of cell labels in cell metadata.
            Could be "all" '''
    )

    parser.add_argument(
        '-cl_col', '--cell_label_column',
        action='store',
        dest='cl_col',
        default='cluster',
        help='Column that contains the labels from "cell_labels"'
    )

    parser.add_argument(
        '-od', '--output_dir',
        action='store',
        dest='od',
        default='mast_mtx',
        help='Basename of output directory.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Load the AnnData file and set matrix to be log data
    adata = sc.read_h5ad(filename=options.h5ad)

    # First subet cells
    if options.cl != 'all':
        cell_label = options.cl.split(',')
        adata = adata[adata.obs[options.cl_col].isin(cell_label)]

    # Check to see if only one condition exists
    if not is_interaction(options.target):
        if adata.obs[options.target].dtype.name == 'category' and len(
                np.unique(adata.obs[options.target].cat.codes)) <= 1:
            raise Exception('There is only 1 condition.')

    # Get normal cp10k counts as well
    adata.layers['cp10k'] = np.expm1(adata.layers['log1p_cp10k'])

    # Run the conversion function.
    for layer in ['counts', 'log1p_cp10k', 'cp10k']:
        adata.X = adata.layers[layer]
        _ = convert_h5ad_R(
            adata,
            out_dir='{}/{}'.format(options.od, layer)
        )

    _ = write_metadata(
        adata,
        adata.obs.columns.tolist(),
        out_dir='{}'.format(options.od)
    )


if __name__ == '__main__':
    main()
