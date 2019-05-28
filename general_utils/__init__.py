
import projection_virtual_side


def sort_asterix(pattern, convert=str):
    """
    Function to create a sorting function for a given pattern using a callable to convert what ever is found with asterix.

    Parameters
    ----------
    pattern : str
        Pattern to sort for.
        If there is no asterix in the pattern, the `convert` will be called on the whole pattern.
        If there is *one* asterix in the pattern, `convert` will be called on what ever is in the asterix.
        If there are multiple asterix in the pattern, `convert` will be call consecutivly (NOT IMPLEMENTED YET).
    convert : Callable, optional
        Some callable which should be applied. Default is `str`.

    Returns
    -------
    func : callable
        Function that converts a given argument.

    TODO
    ----
    * Implement multiple asterix
    * go to index based slicing
    """
    n_asterix = pattern.count('*')
    if n_asterix == 0:
        return lambda x: convert(x)
    elif n_asterix == 1:
        i = pattern.find('*')
        return lambda x: convert(x.replace(pattern[:i], '').replace(pattern[ i +1:], ''))
    else:
        raise UserWarning('Not implemented yet')
        # import re
        # for i in re.finditer('\*', pattern):
