import pandas as pd
import numpy as np


def dataframe_to_mediawiki(dataframe, columns=None, use_index=True):
    """
    Creates a MediaWiki table from a pandas dataframe.
    Parameters
    ----------
    dataframe : pd.DataFrame
    columns : list or None, optional
        Define custom column names.
    use_index : bool, optional
        Creates a column with the index of the pandas table. Default is `True`.

    Returns
    -------
    str_out : str
        MediaWiki table as str
    """
    if columns is None:
        columns = list(dataframe.columns)

    values = dataframe.values.astype(str)
    nrows, ncols = values.shape

    if use_index:
        if dataframe.index.name is None:
            columns = ['index', ] + list(columns)
        else:
            columns = dataframe.index.name + list(columns)
        ncols += 1

    f = np.frompyfunc(lambda x: len(x), 1, 1)
    max_length_col = f(values).astype(int).max(0)
    if use_index:
        max_length_index = max([len(str(i)) for i in df.index])
        max_length_col = [max_length_index, ] + list(max_length_col)
    max_length_col = np.max([max_length_col, f(columns)], axis=0)

    pattern = "\n|-\n"
    pattern += '| {:%ds} ' % max_length_col[0]
    pattern += ('|| {:%ds} ' * (ncols - 1)) % tuple(max_length_col[1:])
    # Output
    str_out = '{| class="wikitable sortable"'
    # HEADER
    str_out += (pattern[:4] + pattern[4:].replace('|', '!')).format(*columns)
    if use_index:
        for r in range(nrows):
            str_out += pattern.format(dataframe.index[r], *values[r])
    else:
        for r in range(nrows):
            str_out += pattern.format(*values[r])

    str_out += '\n|}'
    return str_out