#!/usr/bin/env python


import argparse
import numpy as np
import pandas as pd
import plotnine as plt9
import matplotlib.pyplot as plt
import seaborn as sns

# avoid tk.Tk issues
import matplotlib
matplotlib.use('Agg')
# matplotlib.style.use('ggplot')


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def plot_unity(xdata, ydata, **kwargs):
    mn = min(xdata.min(), ydata.min())
    mx = max(xdata.max(), ydata.max())
    points = np.linspace(mn, mx, 100)
    plt.gca().plot(
        points,
        points,
        color='k',
        marker=None,
        linestyle='--',
        linewidth=1.0
    )


def calculate_expected_pval(df):
    n = len(df)  # retrive the total number of statistical tests made
    increment = 1.0/np.float64(n)  # calculate increment by which sorted
                                   # p vals should increase under the null
    current_expected = np.float64(1.0)  # starting value for 'expected' pvals
    expected = []
    # retrive list of p values in descending order (least significant first)
    df = df.sort_values(by='pvalue', ascending=False)
    for index, row in df.iterrows():
        expected.append(-np.log10(current_expected))
        current_expected -= increment
    df['expected_pval_neglog10'] = expected
    return(df)


def plot_qq(
    df,
    color_var,
    facet_var=None,
    title=''
):
    """
    Inspired by https://www.cureffi.org/2012/08/15/qq-plots-with-matplotlib/
    """
    # retrive pmin, the most significant (i.e. min) p value (for defining
    # the axes)
    axis_max = max(df['pval_neglog10'])

    if facet_var is None:
        pvals = df.groupby(by=color_var).apply(
            calculate_expected_pval
        ).reset_index(level=color_var, drop=True)
    else:
        pvals = df.groupby(by=[color_var, facet_var]).apply(
            calculate_expected_pval
        ).reset_index(level=[color_var, facet_var], drop=True)

    # now plot these two arrays against each other
    n_colors = pvals[color_var].nunique()
    qqplot = plt9.ggplot(pvals, plt9.aes(
        x='expected_pval_neglog10',
        y='pval_neglog10',
        color=color_var
    ))
    qqplot = qqplot + plt9.geom_point(size=0.1, alpha=0.25)
    qqplot = qqplot + plt9.geom_abline(
        slope=1,
        intercept=0,
        color='black',
        linetype='dashed'
    )
    qqplot = qqplot + plt9.theme_bw()
    if n_colors < 9:
        qqplot = qqplot + plt9.scale_colour_brewer(
            palette='Dark2',
            type='qual'
        )
    qqplot = qqplot + plt9.labs(
        x='Expected (-log10 p-value)',
        y='Observed (-log10 p-value)',
        title=title,
        color=''
    )
    qqplot = qqplot + plt9.lims(
        x=(0, axis_max),
        y=(0, axis_max)
    )
    if facet_var is not None:
        qqplot = qqplot + plt9.facet_wrap(
            '~ {}'.format(facet_var),
            ncol=5
        )
    qqplot = qqplot + plt9.theme(
        strip_text=plt9.element_text(size=5),
        axis_text_x=plt9.element_text(angle=-45, hjust=0)
    )
    # set guide legend alpha to 1
    qqplot = qqplot + plt9.guides(
        color=plt9.guide_legend(override_aes={'size': 2.0, 'alpha': 1.0})
    )
    return(qqplot)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Compare results from differential expression tests.
            """
    )

    parser.add_argument(
        '-df', '--dataframe',
        action='store',
        dest='df',
        required=True,
        help='Results dataframe. Required columns: \
            annot_id, pvalue, annot_coef, \
            cell_label, coef_value.'
    )

    parser.add_argument(
        '-compare_columns', '--columns_to_compare',
        action='store',
        dest='columns_to_compare',
        default='de_method,formula_passed,signed_ranking',
        help='Comma separated list of columns to use to aggregate results to \
            compare.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='output_file',
        default='gsea_results_comparison',
        help='Basename of output file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get the parameters
    file = options.df
    cols_compare = options.columns_to_compare.split(',')
    output_file = options.output_file

    df = pd.read_csv(file, sep='\t', low_memory=False)
    df = df.loc[~pd.isna(df['pvalue'])] # filter out unconverged gene sets

    # Add a key based on columns to compare (e.g., differential expression
    # method, condition, covariates)
    # If nan in cols to compare, set to none
    for i in cols_compare:
        df[cols_compare] = df[cols_compare].fillna('none')
    facet_column = 'facet_key'
    df[facet_column] = df[cols_compare].apply(
        lambda row: '\n'.join(row.values.astype(str)),
        axis=1
    )

    # Add neglog10 value
    small_value = np.empty(len(df), dtype=np.float64)
    small_value.fill(0)
    filt = df['pvalue'] == 0.0
    small_value[filt] = np.nanmin(df['pvalue'][np.invert(filt)])  # ** 1.5
    df['pval_neglog10'] = np.log10(df['pvalue'] + small_value) * -1
    df['pval_signedneglog10'] = df['pval_neglog10'] * np.sign(df['annot_coef'])

    # For each combination of columns...
    # 1. Plot p-value distribution
    # 2. Plot qq plot
    # 3. Plot -log10 pvalue
    # 4. Plot signed -log10 pvalue

    # Plot p-value distribution
    gplt_pvd = plt9.ggplot(df, plt9.aes(
        x='pvalue'
    ))
    gplt_pvd = gplt_pvd + plt9.geom_histogram()
    gplt_pvd = gplt_pvd + plt9.theme_bw()
    gplt_pvd = gplt_pvd + plt9.labs(
        x='p-value',
        title=''
    )
    gplt = gplt_pvd + plt9.facet_wrap(
        '~ {}'.format(facet_column),
        ncol=5
    )
    gplt = gplt + plt9.theme(
        strip_text=plt9.element_text(size=5),
        axis_text_x=plt9.element_text(angle=-45, hjust=0)
    )
    n_facets = df[facet_column].nunique()
    gplt.save(
        '{}-pvalue_dist.png'.format(output_file),
        #dpi=300,
        width=6*(n_facets/4),
        height=4*(n_facets/4),
        limitsize=False
    )
    # Facet by celltype label
    gplt = gplt_pvd + plt9.facet_grid(
        '{} ~ {}'.format('cell_label', facet_column)
    )
    gplt = gplt + plt9.theme(
        strip_text=plt9.element_text(size=5),
        axis_text_x=plt9.element_text(angle=-45, hjust=0)
    )
    n_facets_x = df[facet_column].nunique()
    n_facets_y = df['cell_label'].nunique()
    gplt.save(
        '{}-pvalue_dist-split_cell_label.png'.format(output_file),
        #dpi=300,
        width=6*(n_facets_x/4),
        height=4*(n_facets_y/4),
        limitsize=False
    )

    # QQ Plots
    # Sample specific qq-plots
    for group_name, df_group in df.groupby('cell_label'):
        qq_plot = plot_qq(
            df_group,
            facet_column,
            None,
            title='cell label: {}'.format(group_name)
        )
        qq_plot.save(
            '{}-qq-cell_label__{}.png'.format(output_file, group_name),
            #dpi=300,
            width=6,
            height=4,
            limitsize=False
        )
    # Single qq
    qq_plot = plot_qq(
        df,
        facet_column,
        None
    )
    qq_plot.save(
        '{}-qq.png'.format(output_file),
        #dpi=300,
        width=6,
        height=4,
        limitsize=False
    )
    # Single qq faceted by cell_label
    qq_plot = plot_qq(
        df,
        facet_column,
        'cell_label'
    )
    qq_plot.save(
        '{}-qq-split_cell_label.png'.format(output_file),
        #dpi=300,
        width=10,
        height=8,
        limitsize=False
    )

    # Prep to make dataframe to compare p-values across models
    df['row_index'] = df[
        ['annot_id', 'cell_label', 'coef_value']
    ].apply(
        lambda row: '::'.join(row.values.astype(str)),
        axis=1
    )

    # Make a dataframe to compare p-values across models
    df_plt = df.pivot(
        index='row_index',
        columns='facet_key',
        values='pval_neglog10'
    )
    # Compare -log10 p-value
    sns_plt = sns.pairplot(
        df_plt,
        markers='+',
        kind='reg',
        plot_kws={
            'line_kws': {
                'color': '#ffb90f'
            },
            'scatter_kws': {
                'alpha': 0.1,
                'edgecolor': 'b'
            }
        }
    )
    plt.title('-log10(p-value)')
    sns_plt.map_offdiag(plot_unity)
    sns_plt.savefig(
        '{}-pvalue_neglog10.png'.format(output_file)
    )
    plt.close('all')
    df_plt = df_plt.reset_index()
    df_plt[['annot_id', 'cell_label', 'coef_value']] = df_plt[
        'row_index'
    ].str.split('::', expand=True)
    del df_plt['annot_id']
    df_plt = df_plt.groupby('cell_label')
    for group_name, df_group in df_plt:
        del df_group['cell_label']
        sns_plt = sns.pairplot(
            df_group,
            markers='+',
            kind='reg',
            plot_kws={
                'line_kws': {
                    'color': '#ffb90f'
                },
                'scatter_kws': {
                    'alpha': 0.1,
                    'edgecolor': 'b'
                }
            }
        )
        plt.title('-log10(p-value); cell label: {}'.format(group_name))
        sns_plt.map_offdiag(plot_unity)
        sns_plt.savefig(
            '{}-pvalue_neglog10-cell_label__{}.png'.format(
                output_file,
                group_name
            )
        )
        plt.close('all')

    # Make a dataframe to compare signed p-values across models
    df_plt = df.pivot(
        index='row_index',
        columns='facet_key',
        values='pval_signedneglog10'
    )
    # Compare -log10 p-value
    sns_plt = sns.pairplot(
        df_plt,
        markers='+',
        kind='reg',
        plot_kws={
            'line_kws': {
                'color': '#ffb90f'
            },
            'scatter_kws': {
                'alpha': 0.1,
                'edgecolor': 'b'
            }
        }
    )
    plt.title('Signed -log10(p-value)')
    sns_plt.map_offdiag(plot_unity)
    sns_plt.savefig(
        '{}-pvalue_signedneglog10.png'.format(output_file)
    )
    plt.close('all')
    df_plt = df_plt.reset_index()
    df_plt[['annot_id', 'cell_label', 'coef_value']] = df_plt[
        'row_index'
    ].str.split('::', expand=True)
    del df_plt['annot_id']
    df_plt = df_plt.groupby('cell_label')
    for group_name, df_group in df_plt:
        del df_group['cell_label']
        sns_plt = sns.pairplot(
            df_group,
            markers='+',
            kind='reg',
            plot_kws={
                'line_kws': {
                    'color': '#ffb90f'
                },
                'scatter_kws': {
                    'alpha': 0.1,
                    'edgecolor': 'b'
                }
            }
        )
        plt.title('Signed -log10(p-value); cell label: {}'.format(group_name))
        sns_plt.map_offdiag(plot_unity)
        sns_plt.savefig(
            '{}-pvalue_signedneglog10-cell_label__{}.png'.format(
                output_file,
                group_name
            )
        )
        plt.close('all')


if __name__ == '__main__':
    main()
