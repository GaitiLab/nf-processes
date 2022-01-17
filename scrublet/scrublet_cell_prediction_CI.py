import scrublet as scr
import scipy.io
import pandas as pd
import numpy as np
import argparse
import warnings

"""
Read a sparse scRNA gene matrix and output the scrublet prediction for each cell
"""


def main():
    parser = argparse.ArgumentParser(description='Read a sparse gene matrix and predict doublet frequency using '
                                                 'scrublet')
    parser.add_argument('--input_matrix', '-i', type=str, help='Input sparse matrix for scrublet', required=True)
    parser.add_argument('--output_csv', '-o', type=str, help='Output CSV file containing the cell index,'
                                                             'scrublet score, and scrublet prediction', required=True)
    parser.add_argument('--expected_rate', '-e', type=float, help='Expected doublet rate for scrublet', default=0.05)
    parser.add_argument('--min_counts', '-mu', type=float, help='Minimum UMI counts for scrublet', default=2)
    parser.add_argument('--min_cells', '-mc', type=float, help='Minimum cell number for scrublet', default=3)
    parser.add_argument('--gene_variability', '-gv', type=float, help='Minimum gene variability for scrublet',
                        default=85)
    parser.add_argument('--princ_components', '-pc', type=float, help='Minimum number of principal components '
                                                                      'for scrublet', default=30)
    parser.add_argument('--transpose', '-t', action="store_true", help='Transpose the input matrix')

    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

    args = parser.parse_args()

    counts_matrix = scipy.io.mmread(args.input_matrix).T.tocsc()

    input_matrix = counts_matrix.transpose() if args.transpose else counts_matrix

    scrub = scr.Scrublet(input_matrix, expected_doublet_rate=args.expected_rate)

    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=args.min_counts,
                                                              min_cells=args.min_cells,
                                                              min_gene_variability_pctl=args.gene_variability,
                                                              n_prin_comps=args.princ_components)

    sparse_to_df = pd.DataFrame.sparse.from_spmatrix(input_matrix)
    sparse_to_df["cell_index"] = np.arange(sparse_to_df.shape[0]) + 1
    sparse_to_df["doublet_score"] = doublet_scores.tolist()
    sparse_to_df["doublet_prediction"] = predicted_doublets.tolist()

    sparse_to_df[['cell_index', 'doublet_score', 'doublet_prediction']].to_csv(args.output_csv, index=False)


if __name__ == '__main__':
    main()
