import numpy as np

class Scheeme:

    def __init__(self, function: str, time_partition : int, space_partition : int):
        
        assert time_partition == space_partition, \
            f"Expected equal sizes of the space partition and the time partition, got \
                {time_partition}, {space_partition} "

        self._tdim = time_partition
        self._zdim = space_partition
        self._u = np.zeros(shape=(self._tdim, self._zdim))
        self._A = np.zeros(shape=(self._tdim, self._zdim))
        self._d = 0

    def build_sys(self):
        
        # TODO 
        # increase definition of A to encompase A(-1/2*d)^-1A(1/2*d)
        # there should also be a matrix for the coefficients a_ij
        self._A = np.zeros(shape=(self._tdim, self._zdim))
        self._A[np.eye(len(self._A), k=0, dtype='bool')] = 1 - 2*self._d
        self._A[np.eye(len(self._A), k=1, dtype='bool')] = self._d
        self._A[np.eye(len(self._A), k=-1, dtype='bool')] = self._d
        return 0

    def calculate(self, S, r, sigma, K, T):

        for t in range(self._tdim):
            
            if t == 0:
                continue

            # TODO 
            # provide calculation for aquiring Z
            Z = ...
            self._d = T * self._tdim**2/ (Z**2 * self._zdim)
            self.build_sys()
            self._u[1:,t+1] = self._A*self._u[1:,t]




