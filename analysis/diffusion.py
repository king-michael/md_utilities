import numpy as np
from scipy.optimize import curve_fit
from collections import deque
import MDAnalysis

try:
    from progress_reporter import ProgressReporter_
except:
    from ..fallbacks import ProgressReporter_


def calculate_msd(ag, n_min=1, n_max_number_steps=np.inf, dt=1000):
    """
    Parameter
    ---------
    ag :
        atomgroup
    n_min : int
        Minimum number of intervals to contribute to diffusivity.
    n_max_number_steps : int
        Maximum number of intervals to contribute to diffusivity.
        Spans together with `n_min` the window which is moved along the trajectory.
        Use `np.inf` to use the maximum possible length (`int(np.floor(n_frames-1) / 2)`).
        Default is `np.inf`.
    dt: float
        Time step between frames. [fs]

    Returns
    -------
    time_vec : np.darray
        Time vector
    xmsd : np.ndarray
        MSD vector
    xmsd_var : np.ndarray
        Variance of the MSD

    Examples
    --------
    >>> u = MDAnalysis.Universe('water_30A_mW.psf', 'trajectory.1.dcd')
    >>> ag = u.atoms.select_atoms('type mW')
    >>> dt = 5000 # time in fs
    >>> time_vec, xmsd, xmsd_var = calculate_msd(ag, n_min=1, n_max_number_steps=np.inf, dt=dt)
    """
    assert n_min >= 1, 'Need at least one frame!'

    n_frames = ag.universe.trajectory.n_frames
    n_atoms = ag.n_atoms

    n_origin = int(np.floor(n_frames - 1) / 2)
    #  maximum number of intervals to contribute to diffusivity
    n_max = min(n_max_number_steps, n_origin)
    n_origin = n_frames - n_max
    n_min = min(n_min, n_max)

    n_window = n_max - n_min
    #  store mean square displacements in xmsd
    time_vec = np.arange(dt * (n_min + 1), dt * (n_max + 1), dt, dtype=np.float64)
    xmsd = np.zeros((n_window, 3), dtype=np.float64)
    xmsd2 = np.zeros((n_window, 3), dtype=np.float64)

    # temporary storage
    deque_coords = deque(maxlen=n_max)
    pgr = ProgressReporter_()
    pgr.register(n_frames - n_window)
    for i, ts in enumerate(ag.universe.trajectory):
        deque_coords.append(ts.positions[ag._ix])
        if i >= (n_max - 1):
            pgr.update(1)
            coords_t0 = deque_coords.popleft()
            msd_current = np.sum(np.power(np.asarray(deque_coords)[n_min - 1:] - coords_t0, 2), axis=1)
            xmsd += msd_current
            xmsd2 += np.power(msd_current, 2)
    pgr.finish()
    xmsd /= n_atoms * n_origin
    xmsd2 /= n_atoms * n_origin

    xmsd_var = xmsd2 - np.power(xmsd, 2)

    return time_vec, xmsd, xmsd_var


def calculate_diffusion_coefficent(time, msd, return_per_dim=False, return_fit=False, verbose=False):
    """
    Use the follwing units: Time [fs] and MSD : [A^2]

    Parameter
    ---------
    time : np.ndarray
        Time [fs] 1-D time vector
    msd : np.ndarray
        MSD [A^2] array. Can be 1-D or 2-D with shape `(n_points,)` or `(n_points, n_dim)`.
    return_per_dim : bool, optional
        Returns also the diffusion coefficent in the different directions.
        Default is `False`.
    return_fit : bool, optional
        Returns also the parameter and covariance matrix of the linear fit.
        Default is `False`.
    verbose : bool, optional
        Turns on verbosity. Provides informations about the fit and the resulting diffusion coefficents.
        Default is `False`.

    Returns
    -------
    (Dav, Davsd) : tuple[float, float]
        Average diffusion coefficent and standard deviation of it.
    (D, Dsd) : tuple[np.darray, np.darray], optional
        Diffusion coefficent in each dimension and standard deviation of them.
    (popt, pcov) : tuple[np.darray, np.darray], optional
        Parameter and covariance matrix of the linear fit.

    Examples
    --------
    Average diffusion coefficient and standard deviation
    >>> (Dav, Davsd) = calculate_diffusion_coefficent(time, msd, verbose=True)

    Average diffusion coefficient and standard deviation for sliced data.
    >>> xslice = slice(10, None, None)
    >>>> (Dav, Davsd) = calculate_diffusion_coefficent(time[xslice], msd[xslice], verbose=True)

    Average diffusion coefficient and the diffiusion coefficent per dimension
    >>>> (Dav, Davsd), (D, Dsd) = calculate_diffusion_coefficent(time, msd, return_per_dim=True)

    Average diffusion coefficient and fitting parameter
    >>>> (Dav, Davsd), (popt, pcov) = calculate_diffusion_coefficent(time, msd, return_fit=True)
    """
    # convert_units = 1.0e-5, #convert A^2/fs to m^2/sec
    # convert_units = 0.1,    #convert A^2/fs to cm^2/sec
    convert_units = 0.1
    linear = lambda x, m, n: m * x + n

    if msd.ndim == 1:
        msd = np.atleast_2d(msd).T

    N, ndims = msd.shape

    if verbose:
        print('Data')
        print('  n_points : {}'.format(N))
        print('  n_dims   : {}'.format(ndims))
    assert N == time.size, \
        "MSD do not match time.\n time : {}\nMSD : {}".format(time.shape, msd.shape)

    popt = np.empty((ndims, 2), dtype=np.float64)
    pcov = np.empty((ndims, 2, 2), dtype=np.float64)

    for i in range(ndims):
        popt[i], pcov[i] = curve_fit(linear, time, msd[:, i])

    if verbose:
        print('\nFitting results')
        for d, p, c in zip('xyz', popt, pcov):
            print('  {:s} : slope     ={: e} (\u00B1 {:e})\n'.format(d, p[0], np.sqrt(c[0, 0])) +
                  '      intercept ={: e} (\u00B1 {:e})\n'.format(p[1], np.sqrt(c[1, 1])))

    D = 0.5 * popt[:, 0] * convert_units
    Dsd = np.sqrt(pcov[:, 0, 0])

    Dav = np.mean(D)
    Davsd = np.std(D)

    if verbose:
        print('Diffusion Coefficent')
        for d, d_dim, ssd_dim in zip('xyz', D, Dsd):
            print('  {:3s} : {: e} (\u00B1 {:e}) cm^2/sec'.format(d, d_dim, ssd_dim))
        print('  {:3s} : {: e} (\u00B1 {:e}) cm^2/sec'.format('avg', Dav, Davsd))

    if return_fit:
        if return_per_dim:
            return (Dav, Davsd), (D, Dsd), (popt, pcov)
        else:
            return (Dav, Davsd), (popt, pcov)
    if return_per_dim:
        return (Dav, Davsd), (D, Dsd)

    return Dav, Davsd