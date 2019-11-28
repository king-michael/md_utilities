import numpy as np


class ReaderPolymorphic(object):
    def __init__(self, filename):
        self.filename = filename
        self.read(filename)

    def read(self, filename):
        with open(filename, 'r') as fp:
            while True:
                line = fp.readline().strip()
                if line[0] == '#':
                    continue
                else:
                    break
            # types / self.eta
            self.ntypes, self.eta = map(int, line.split())  # pair_polymorphic.cpp, line 586 / 591
            self.npairs = self.ntypes * (self.ntypes + 1) // 2  # pair_polymorphic.cpp, line 637
            self.ntriple = self.ntypes * self.ntypes * self.ntypes  # pair_polymorphic.cpp, line 638
            assert self.eta == 0, "self.eta has to be 0 for sw"
            # tpye section
            self.types = []
            self.masses = []
            self.atomic_number = []
            for i in range(self.ntypes):  # pair_polymorphic.cpp, line 647
                line = fp.readline().strip()
                self.atomic_number.append(line.split()[0])
                self.masses.append(line.split()[1])
                self.types.append(line.split()[2])

            self.pairs = [(t, t) for t in self.types]
            self.pairs.extend([(t1, t2) for i, t1 in enumerate(self.types)
                               for j, t2 in enumerate(self.types)
                               if i < j])
            self.triples = [(t1, t2, t3) for t1 in self.types
                            for t2 in self.types
                            for t3 in self.types]

            # bins
            self.nr, self.ng, self.nx, self.maxX = map(int, map(float,
                                                                fp.readline().split()))  # pair_polymorphic.cpp, line 623 - 631
            # Cutoff section
            self.cutoffs = []
            self.xi = []
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 647
                cut, xi = map(float, fp.readline().split())
                self.cutoffs.append(cut)
                self.xi.append(xi)
            self.cutmax = 0
            for i in range(self.npairs):
                cut = self.cutoffs[i]
                if cut > self.cutmax: self.cutmax = cut

            lines = fp.read()  # read the rest
            lines = [float(i) for i in lines.split()]

            offset = 0
            # U
            tmp = []
            range_U = []
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 671
                tmp.append(lines[offset:offset + self.nr])  # pair_polymorphic.cpp, line 675
                offset += self.nr
                range_U.append([0, self.cutoffs[i]])
            self.U = np.array(tmp)

            # V
            tmp = []
            range_V = []
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 681
                tmp.append(lines[offset:offset + self.nr])  # pair_polymorphic.cpp, line 684
                offset += self.nr
                range_V.append([0, self.cutoffs[i]])
            self.V = np.array(tmp)
            # W
            tmp = []
            range_W = []
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 690
                tmp.append(lines[offset:offset + self.nr])  # pair_polymorphic.cpp, line 693
                offset += self.nr
                range_W.append([0, self.cutoffs[i]])
            self.W = np.array(tmp)
            # P
            tmp = []
            range_P = []
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 699
                tmp.append(lines[offset:offset + self.nr])  # pair_polymorphic.cpp, line 702
                offset += self.nr
                range_P.append([-self.cutmax, self.cutmax])
            self.P = np.array(tmp)
            # G
            tmp = []
            range_G = []
            for i in range(self.ntriple):  # pair_polymorphic.cpp, line 710
                tmp.append(lines[offset:offset + self.ng])  # pair_polymorphic.cpp, line 713
                offset += self.ng
                range_G.append([-1, 1])
            self.G = np.array(tmp)
            # F
            tmp = []
            range_F = []
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 721
                tmp.append(lines[offset:offset + self.nx])  # pair_polymorphic.cpp, line 724
                offset += self.nx
                range_F.append([0, self.maxX])
            self.F = np.array(tmp)


class WriterPolymorphic():
    values_per_line = 5
    verbose = True

    def __init__(self, list_types, list_cutoff, num_bins):
        self.list_types = list_types
        self.ntypes = len(list_types)
        self.types = [t[2] for t in list_types]
        self.pairs = [(t, t) for t in self.types]
        self.pairs.extend([(t1, t2) for i, t1 in enumerate(self.types)
                           for j, t2 in enumerate(self.types)
                           if i < j])
        self.triples = [(t1, t2, t3) for t1 in self.types
                        for t2 in self.types
                        for t3 in self.types]

        self.list_cutoff = list_cutoff
        self.num_bins = num_bins

    def check_parameter(self):
        n_pairs = self.ntypes * (self.ntypes + 1) // 2
        n_triples = self.ntypes * self.ntypes * self.ntypes

        assert len(self.list_types) == self.ntypes, "Number of types and self.list_types don't match"
        assert len(self.list_cutoff) == n_pairs, "Number of pairs and self.list_cutoff don't match"

        assert len(self.pairs) == n_pairs, "Number of pairs and self.pairs don't match"
        assert len(self.triples) == n_triples, "Number of triples and self.triples don't match"

        assert hasattr(self, 'potential_U'), \
            "Define self.potential_U with shape ({}, {})".format(n_pairs, self.num_bins[0])
        assert hasattr(self, 'potential_V'), \
            "Define self.potential_V with shape ({}, {})".format(n_pairs, self.num_bins[0])
        assert hasattr(self, 'potential_W'), \
            "Define self.potential_W with shape ({}, {})".format(n_pairs, self.num_bins[0])
        assert hasattr(self, 'potential_P'), \
            "Define self.potential_P with shape ({}, {})".format(n_pairs, self.num_bins[0])
        assert hasattr(self, 'potential_G'), \
            "Define self.potential_G with shape ({}, {})".format(n_triples, self.num_bins[1])
        assert hasattr(self, 'potential_F'), \
            "Define self.potential_F with shape ({}, {})".format(n_pairs, self.num_bins[2])

    def write(self, fname):
        self.check_parameter()
        with open(fname, "w") as fp:
            self._write_header(fp)
            self._write_body(fp)

    def _write_header(self, fp):

        eta = 0

        n_pairs = self.ntypes * (self.ntypes + 1) // 2

        # Write Comment
        STR = "### created from Michael King at {}\n".format(time.strftime("%Y-%m-%d"))
        # write ntypes, eta
        STR += " {} {}\n".format(self.ntypes, eta)
        # write atomic_number atomic-mass element (ntypes)
        for i in range(self.ntypes):
            entry = self.list_types[i]
            STR += " {} {} {}\n".format(*entry)
        # write number of bins
        STR += "{} {} {} {:e}\n".format(*num_bins)
        # wrtie cut xi (ntypes*(ntypes+1)/2)
        for i in range(n_pairs):
            xi = 0
            STR += " {} {}\n".format(self.list_cutoff[i], xi)
        fp.write(STR)
        if self.verbose:
            print(STR)

    def _write_potential(self, fp, pot):
        "writes a potential"
        fmt_line = "{: .8e} " * self.values_per_line
        fmt_line = fmt_line[:-1] + '\n'
        pot_re = pot.reshape((-1, self.values_per_line))
        for i, line in enumerate(pot_re):
            fp.write(fmt_line.format(*line))

    def calculate_VW(self, X, eps, sig, a, lam, gam):
        """Stillinger-Weber (distance-based splitted)
        Y=sqrt(lam*eps) * exp( gam*a / (r - a * sig) )
        return: Y
        """
        X = np.array(X)
        Y = np.zeros((len(X)))
        rcut = sig * a
        prefactor = np.sqrt(lam * eps)
        if prefactor == 0.0:
            return Y
        expfactor = gam * sig
        tmp = np.where(X < rcut)
        Y[tmp] = prefactor * np.exp(expfactor / (X[tmp] - rcut))
        return Y

    def calculate_U(self, X, eps, sig, a, A, B, p, q):
        """Stillinger-Weber (distance-based splitted)
        Y=sqrt(lam*eps) * exp( gam*a / (r - a * sig) )
    return: Y"""
        X = np.array(X)
        Y = np.zeros((len(X)))
        rcut = sig * a
        prefactor = A * eps
        if prefactor == 0.0:
            return Y
        tmp = np.where(X < rcut)
        Y[tmp] = prefactor * (B * np.power(sig / X[tmp], p) - np.power(sig / X[tmp], q)) * np.exp(sig / (X[tmp] - rcut))
        return Y

    def Chebyshev(self, m, c):
        """Chebyshev polynomials
        arcos(cos(n\theta))=Chebyshev(m,cos(theta))
        tn   = 1.0
        tn_1 = 1.0
        tn_2 = 0.0
        tn_2 = c
        for i in np.arange(1,m+1):
            tn = 2*c*tn_1 - tn_2;
            tn_2 = tn_1;
            tn_1 = tn;
        return tn
        """
        tn = 1.0
        tn_1 = 1.0
        tn_2 = 0.0
        tn_2 = c
        for i in np.arange(1, m + 1):
            tn = 2 * c * tn_1 - tn_2;
            tn_2 = tn_1;
            tn_1 = tn;
        return tn

    # def calculate_G(self, C, B, n, nbins):
    #     # same as  x= np.cos(theta) ; y =  C*(1-B*(-1)**n*np.cos(n*theta))
    #     B = int(B)
    #     n = int(n)
    #     x = np.linspace(-1, 1, nbins)
    #     tn = self.Chebyshev(n, x)
    #     Y = C * (1 - B * np.power(-1, n) * tn)
    #     return Y

    def calculate_G(self, costheta, nbins):
        "Stillinger Weber G(theta) term"
        t = np.linspace(-1, 1, nbins)
        return np.power(t - costheta, 2)


    def calculate_potentials(self, params_SW, params_G):
        n_pairs = self.ntypes * (self.ntypes + 1) // 2
        n_triples = self.ntypes * self.ntypes * self.ntypes

        assert params_SW.shape[0] == n_pairs, "Number of pairs and params_SW don't match"
        assert params_G.shape[0] == n_triples, "Number of pairs and params_G don't match"

        assert len(set([j for i in params_SW.index for j in i])) == self.ntypes, \
            'There are duplciates in params_SW'
        assert len(set([j for i in params_G.index for j in i])) == self.ntypes, \
            'There are duplciates in params_G'
        index_list = params_SW.index

        def fix_inf(pot):
            pot[pot == np.inf] = pot[pot != np.inf].max()
            return pot

        self.potential_U = np.zeros((n_pairs, self.num_bins[0]))
        self.potential_V = np.zeros((n_pairs, self.num_bins[0]))
        self.potential_W = np.zeros((n_pairs, self.num_bins[0]))
        self.potential_P = np.zeros((n_pairs, self.num_bins[0]))
        self.potential_F = np.zeros((n_pairs, self.num_bins[2]))

        for i, p in enumerate(self.pairs):
            x = np.linspace(0, self.list_cutoff[i], self.num_bins[0])

            pot = self.calculate_U(x, *params_SW.loc[p][['epsilon', 'sigma', 'a', 'A', 'B', 'p', 'q']])
            self.potential_U[i] = fix_inf(pot)

            pot = self.calculate_VW(x, *params_SW.loc[p][['epsilon', 'sigma', 'a', 'lambda', 'gamma']])
            self.potential_V[i] = fix_inf(pot)

            pot = self.calculate_VW(x, *params_SW.loc[p][['epsilon', 'sigma', 'a', 'lambda', 'gamma']])
            self.potential_W[i] = fix_inf(pot)

            self.potential_P[i] = np.ones(self.num_bins[0])

            self.potential_F[i] = np.linspace(-1.00000000e-20, -self.num_bins[3], self.num_bins[2])

        self.potential_G = np.zeros((n_triples, self.num_bins[1]))
        for i, t in enumerate(self.triples):
            pot = self.calculate_G(*params_G.loc[t], nbins=self.num_bins[1])
            self.potential_G[i] = pot

    def _write_body(self, fp):
        n_pairs = self.ntypes * (self.ntypes + 1) // 2
        n_triples = self.ntypes * self.ntypes * self.ntypes

        for i in range(n_pairs):  # U
            self._write_potential(fp, self.potential_U[i])
        for i in range(n_pairs):  # V
            self._write_potential(fp, self.potential_V[i])
        for i in range(n_pairs):  # W
            self._write_potential(fp, self.potential_W[i])
        for i in range(n_pairs):  # P
            self._write_potential(fp, self.potential_P[i])
        for i in range(n_triples):  # G
            self._write_potential(fp, self.potential_G[i])
        for i in range(n_pairs):  # F
            self._write_potential(fp, self.potential_F[i])