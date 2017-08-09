import pandas as pd
import numpy as np
import bisect


def compute_p_value(scores, null_p_values):
    """Get the p-value for each score by examining the list null distribution
    where scores are obtained by a certain probability.

    NOTE: uses score2pval function

    Parameters
    ----------
    scores : pd.Series
        series of observed scores
    null_p_values: pd.Series
        Empirical null distribution, index are scores and values are p values

    Returns
    -------
    pvals : pd.Series
        Series of p values for scores
    """
    num_scores = len(scores)
    pvals = pd.Series(np.zeros(num_scores))
    null_p_val_scores = list(reversed(null_p_values.index.tolist()))
    #null_p_values = null_p_values.ix[null_p_val_scores].copy()
    null_p_values.sort_values(inplace=True, ascending=False)
    pvals = scores.apply(lambda x: score2pval(x, null_p_val_scores, null_p_values))
    return pvals


def score2pval(score, null_scores, null_pvals):
    """Looks up the P value from the empirical null distribution based on the provided
    score.

    NOTE: null_scores and null_pvals should be sorted in ascending order.

    Parameters
    ----------
    score : float
        score to look up P value for
    null_scores : list
        list of scores that have a non-NA value
    null_pvals : pd.Series
        a series object with the P value for the scores found in null_scores

    Returns
    -------
    pval : float
        P value for requested score
    """
    # find position in simulated null distribution
    pos = bisect.bisect_right(null_scores, score)

    # if the score is beyond any simulated values, then report
    # a p-value of zero
    if pos == null_pvals.size and score > null_scores[-1]:
        return 0
    # condition needed to prevent an error
    # simply get last value, if it equals the last value
    elif pos == null_pvals.size:
        return null_pvals.iloc[pos-1]
    # normal case, just report the corresponding p-val from simulations
    else:
        return null_pvals.iloc[pos]


def cummin(x):
    """A python implementation of the cummin function in R"""
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    """A python implementation of the Benjamani-Hochberg FDR method.
    This code should always give precisely the same answer as using
    p.adjust(pval, method="BH") in R.
    Parameters
    ----------
    pval : list or array
        list/array of p-values
    Returns
    -------
    pval_adj : np.array
        adjusted p-values according the benjamani-hochberg method
    """
    pval_array = np.array(pval)
    sorted_order = np.argsort(pval_array)
    original_order = np.argsort(sorted_order)
    pval_array = pval_array[sorted_order]

    # calculate the needed alpha
    n = float(len(pval))
    pval_adj = np.zeros(int(n))
    i = np.arange(1, int(n)+1, dtype=float)[::-1]  # largest to smallest
    pval_adj = np.minimum(1, cummin(n/i * pval_array[::-1]))[::-1]
    return pval_adj[original_order]
