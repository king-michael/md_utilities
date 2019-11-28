import numpy as np

class PolymorphicFile(object):
    def __init__(self, filename):
        self.filename = filename
        self.read(filename)
    def read(self, filename):
        with open(filename, 'r') as fp:
            while True:
                line = fp.readline().strip()
                if line[0] == '#': continue
                else: break
            # types / self.eta
            self.ntypes , self.eta = map(int,line.split()) # pair_polymorphic.cpp, line 586 / 591
            self.npairs = self.ntypes*(self.ntypes+1)//2 # pair_polymorphic.cpp, line 637
            self.ntriple = self.ntypes * self.ntypes * self.ntypes # pair_polymorphic.cpp, line 638
            assert self.eta == 0, "self.eta has to be 0 for sw"
            # tpye section
            self.list_types=[]
            for i in range(self.ntypes): # pair_polymorphic.cpp, line 647
                line = fp.readline().strip()
                self.list_types.append(())
            # bins
            self.nr,self.ng,self.nx,self.maxX  = map(int,map(float, fp.readline().split()))  # pair_polymorphic.cpp, line 623 - 631
            # Cutoff section
            self.list_cutoff= []
            for i in range(self.npairs): # pair_polymorphic.cpp, line 647
                self.list_cutoff.append([float(i) for i in fp.readline().split()])

            self.cutmax = 0
            for i in range(self.npairs):
                cut = self.list_cutoff[i][0]
                if cut > self.cutmax: self.cutmax = cut

            lines = fp.read() # read the rest
            lines = [float(i) for i in lines.split()]

            offset = 0
            # U
            tmp = []
            range_U=[]
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 671
                tmp.append(lines[offset:offset+self.nr]) # pair_polymorphic.cpp, line 675
                offset+=self.nr
                range_U.append([0,self.list_cutoff[i][0]])
            self.U=np.array(tmp)

            # V
            tmp = []
            range_V=[]
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 681
                tmp.append(lines[offset:offset+self.nr]) # pair_polymorphic.cpp, line 684
                offset+=self.nr
                range_V.append([0,self.list_cutoff[i][0]])
            self.V=np.array(tmp)
            # W
            tmp = []
            range_W=[]
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 690
                tmp.append(lines[offset:offset+self.nr]) # pair_polymorphic.cpp, line 693
                offset+=self.nr
                range_W.append([0,self.list_cutoff[i][0]])
            self.W=np.array(tmp)
            # P
            tmp = []
            range_P=[]
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 699
                tmp.append(lines[offset:offset+self.nr]) # pair_polymorphic.cpp, line 702
                offset+=self.nr
                range_P.append([-self.cutmax,self.cutmax])
            self.P=np.array(tmp)
            # G
            tmp = []
            range_G=[]
            for i in range(self.ntriple):  # pair_polymorphic.cpp, line 710
                tmp.append(lines[offset:offset+self.ng]) # pair_polymorphic.cpp, line 713
                offset+=self.ng
                range_G.append([-1,1])
            self.G=np.array(tmp)
            # F
            tmp = []
            range_F=[]
            for i in range(self.npairs):  # pair_polymorphic.cpp, line 721
                tmp.append(lines[offset:offset+self.nx]) # pair_polymorphic.cpp, line 724
                offset+=self.nx
                range_F.append([0,self.maxX])
            self.F=np.array(tmp)