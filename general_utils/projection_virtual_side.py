"""
Python function to add a virtual side to a trajectory of a molecule.
"""
import numpy as np

def get_projections(A, B, C, D, verbose=False):
    """
    Get the projection values for D, as function of points A,B,C

    D := a*AB + b*CB + t * n

    Parameters
    ----------
    A : numpy.darray
        3D point
    B : numpy.darray
        3D point
    C : numpy.darray
        3D point
    D : numpy.darray
        3D point of group of intrest
    verbose : bool
        Default is False.

    Returns
    -------
    a, b, t
    """
    precision = 9
    AB = B - A  # vector: A -> B
    CB = B - C  # vector: C -> B
    # BEGIN planevector
    n_unnorm = np.cross(AB, CB)  # orthogonal vector to AB, CB
    l_n = np.linalg.norm(n_unnorm)
    n = n_unnorm / l_n  # normal vector
    dist_0 = np.dot(n, A)  # distance to the origin
    #     assert dist_0 == np.dot(n,B)
    #     assert dist_0 == np.dot(n,C)
    # END planevector

    if verbose: print("<x,{}> = {}".format(n, dist_0))
    # BEGIN distance point,plane
    t = (np.dot(D, n) - dist_0) / np.dot(n, n)
    # END distance point,plane
    # BEGIN projection point on plane
    proj_point = D - t * n
    S = proj_point - B
    MATRIX = np.array((AB, CB)).T
    MATRIX22 = []
    STMP = []
    COUNT = 0
    for i, line in enumerate(MATRIX):
        if COUNT < 2:
            if np.round(line[0], precision) == np.round(line[1], precision)\
                    and np.round(line[0], precision) == 0:
                pass
            else:
                MATRIX22.append(line)
                STMP.append(S[i])
                COUNT += 1
    MATRIX22 = np.array(MATRIX22)
    MATRIX22_inv = np.linalg.inv(MATRIX22)
    a, b = np.dot(MATRIX22_inv, STMP)
    # END projection point on plane
    if verbose:
        print("D := {}*AB + {}*CB + {}*n".format(a, b, t))
        print("AB: {}".format(AB))
        print("CB: {}".format(CB))
        print("n:  {}".format(n))
    return (a, b, t)


def project(A, B, C, PROJECTIONVALUES):
    """
    Function to calculate the projected coordinates in dependecy of `A, B, C` and `PROJECTIONVALUES`
    PROJECTIONVALUES = (a, b, t)

    D := a*AB + b*CB + t * n

    Parameters
    ----------
    A : numpy.darray
        3D point
    B : numpy.darray
        3D point
    C : numpy.darray
        3D point
    PROJECTIONVALUES
        (a, b, t) returned by `get_projections(A, B, C, D)`

    Returns
    -------
    D : numpy.darray
        3D point
    """
    a, b, t = PROJECTIONVALUES
    AB = B - A  # vector: A -> B
    CB = B - C  # vector: C -> B
    n_unnorm = np.cross(AB, CB)  # orthogonal vector to AB, CB
    l_n = np.linalg.norm(n_unnorm)
    n = n_unnorm / l_n  # normal vector
    dist_0 = np.dot(n, A)  # distance to the origin
    assert dist_0 == np.dot(n, B)
    assert dist_0 == np.dot(n, C)
    D = B + a * AB + b * CB + (t) * n
    return D

if __name__ == '__main__':

    # testcase
    A = np.array([0, 0, 0])
    B = np.array([0, 0, 1])
    C = np.array([0, 1, 0])
    D = np.array([2, 0, 0])
    # B  C
    # | /
    # A ----D
    #
    PROJECTIONVALUES = get_projections(A, B, C, D)
    print("Projection values: {}".format(PROJECTIONVALUES))

    D_test = project(A, B, C, PROJECTIONVALUES)
    np.testing.assert_array_almost_equal(D_test, [2,0,0])
    print("-" * 40 + "\n" + "D should be at [2,0,0]")
    print("D projected: {}".format(D_test))

    # testcase it rotated
    A = np.array([0, 0, 0])
    B = np.array([1, 0, 0])
    C = np.array([0, 1, 0])
    # D should be at [0,0,-2]
    #   C
    #  /
    # A ---B
    # |
    # D
    D_test = project(A, B, C, PROJECTIONVALUES)
    np.testing.assert_array_almost_equal(D_test, [0,0,-2])
    print("-" * 40 + "\n" + "D should be at [0,0,-2]")
    print("D projected: {}".format(D_test))

    # testcase it moved and rotated
    A = np.array([1, 0, 1])
    B = np.array([2, 0, 1])
    C = np.array([1, 1, 1])
    # D should be at [1,0,-1]
    #     C
    #    /
    # --A ---B
    # |  |
    # 0   D
    D_test = project(A, B, C, PROJECTIONVALUES)
    np.testing.assert_array_almost_equal(D_test, [1,0,-1])
    print("-" * 40 + "\n" + "D should be at [1,0,-1]")
    print("D projected: {}".format(D_test))