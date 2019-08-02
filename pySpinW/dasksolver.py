import dask.array as da
import numpy as np
#import cupy
from scipy.linalg import cholesky, eig, LinAlgError


class Hamiltonian:

    def __init__(self, client, hamIn, chunks=1000):
        self.client = client
        self.nHKL = hamIn.shape[0]
        self.nMagExt = hamIn.shape[1]
        self.ham = da.from_array(hamIn, chunks=(chunks, self.nMagExt, self.nMagExt))
        self.storage = {
            'cho': None,
            'V': None,
            'omega': None
        }

    def compute(self):
        return self.client.gather(self.storage['V']), self.client.gather(self.storage['omega'])



    def hermitSolve(self):
        self.storage['cho'] = self.cholesky_gu(self.ham)
        self.storage['V'], self.storage['omega'] = self.hermitian(self.storage['cho'])
        self.storage['V'] = self.client.persist(self.storage['V'])
        self.storage['omega'] = self.client.persist(self.storage['omega'])

    def nonHermitSolve(self):
        raise NotImplementedError

    @staticmethod
    def cholesky_gu(x):
        if x is None:
            raise ValueError
        @da.as_gufunc(signature="(j,j)->(j,j)", output_dtypes=np.dtype('c16'), vectorize=True)
        def f_gu(block):
            try:
                result = cholesky(block, check_finite=False, overwrite_a=False)
            except LinAlgError:
                toll = eig(block)
                toll = np.sort(toll)
                toll = np.abs(toll[0])
                toll = toll*np.sqrt(block.shape[0]*2)*4
                result = cholesky(block + np.eye(block.shape[0])*toll, check_finite=False, overwrite_a=False)
            return result
        return f_gu(x)

    def hermitian(self, x):
        gCommd = np.ones(self.nMagExt)
        gCommd[int(self.nMagExt / 2):] = -1
        gComm = gCommd * np.eye(self.nMagExt)
        @da.as_gufunc("(j,j)->(j,j),(j)", output_dtypes=2*(np.dtype('c16'), ), vectorize=True)
        def eigsolver(block):
            # try:
            temp = block * gComm * block.transpose()
            # except:
            #     return -1
            temp = 0.5 * (temp + temp.transpose())
            omega, V = eig(temp)
            ind = np.argsort(omega.real)
            V = V[:, ind]
            omega = omega[ind]
            # TODO real part of temp should be sorted and also applied to omega
            V = np.linalg.inv(block) * V * np.diag(np.sqrt(gCommd * omega))
            return V, omega
        return eigsolver(x)
