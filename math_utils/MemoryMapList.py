import numpy as np
from itertools import chain


class MemoryMapList(object):
    def __init__(self, npy_files):
        """
        Class to combine multiple MemoryMap files into one.

        Parameter
        ---------
        npy_files : str or List[str] or List[nd.array]
            List of npy files. Can also handle a list of arrays

        Example
        -------
        >>> array = np.vstack([np.arange(20), np.arange(20)]).T
        >>> array_split = np.array_split(array, 4)
        >>> stack = MemoryMapList(array_split)

        test get item
        >>> np.testing.assert_array_equal( np.vstack([stack[i] for i in range(len(stack))]), array )

        test iteration
        >>> np.testing.assert_array_equal(np.vstack([i for i in stack]), array)

        test slice
        >>> np.testing.assert_array_equal(stack[0:len(array)], array)
        >>> np.testing.assert_array_equal(stack[None:None], array)
        >>> np.testing.assert_array_equal(stack[4:16], array[4:16])

        Not implemented
        >>> np.testing.assert_array_equal(stack[None:None:2], array[::2])

        """
        if isinstance(npy_files, str):
            npy_files = [npy_files, ]
        if isinstance(npy_files[0], np.ndarray):
            self._list_arrays = npy_files
        else:
            self._npy_files = npy_files
            self._list_arrays = [np.load(f, mmap_mode='r') for f in npy_files]
        self._shape = np.array([a.shape for a in self._list_arrays], dtype=np.int)
        assert np.unique(self._shape[:, 1:], axis=0).size == 1
        self.shape = (sum(self._shape[:, 0]),) + tuple(self._shape[0, 1:].tolist())

        self._shape_cumsum = np.concatenate([[0], np.cumsum(self._shape[:, 0])])

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, key):
        if isinstance(key, (slice,)):
            start, stop, step = key.start, key.stop, key.step
            if not step is None or step == 1:
                raise NotImplementedError("Not implemented steps")
            if start is None:
                start = 0
            if stop is None:
                stop = self.shape[0]

            istart = np.searchsorted(self._shape_cumsum, start, side='right') - 1
            iend = np.searchsorted(self._shape_cumsum, stop, side='left') - 1

            idiff = (iend - istart)
            # same trajectory
            if idiff == 0:
                offset = self._shape_cumsum[istart]
                key1 = slice(start - offset, stop - offset, step)

                return self._list_arrays[istart][key1]

            # consecutive trajectory
            elif idiff == 1:
                list_slices = [
                    (istart, slice(start - self._shape_cumsum[istart], None, step)),
                    (iend, slice(None, stop - self._shape_cumsum[iend], step)),
                ]

                return np.concatenate([self._list_arrays[i][s] for i, s in list_slices])

            # 3 or more trajectories
            else:
                list_slices = [
                    (istart, slice(start - self._shape_cumsum[istart], None, step)),
                ]
                for i in range(istart + 1, iend):
                    list_slices.append((i, slice(None, None, step)))
                list_slices.append((iend, slice(None, stop - self._shape_cumsum[iend], step)))

                return np.concatenate([self._list_arrays[i][s] for i, s in list_slices])

        # single item
        f = np.searchsorted(self._shape_cumsum, key, side='right') - 1
        i = key - self._shape_cumsum[f]

        return self._list_arrays[f][i]

    def __iter__(self):
        return chain(*[iter(a) for a in self._list_arrays])
