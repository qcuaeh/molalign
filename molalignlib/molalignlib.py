import numpy as np
from ase import Atoms

try:
    from .f2py_molalignlib import molalignlib
except ModuleNotFoundError:
    from os import path
    from subprocess import Popen, PIPE, STDOUT
    moduledir = path.dirname(path.abspath(__file__))
    with Popen(['./compile.sh', '-all', '-pic', '.', '.'], cwd=moduledir, stdout=PIPE, stderr=STDOUT, bufsize=0) as p:
        for line in p.stdout:
            print(line.decode('utf-8').rstrip())
    with Popen(['./link_pythonlib.sh', '.', 'f2py_molalignlib'], cwd=moduledir, stdout=PIPE, stderr=STDOUT, bufsize=0) as p:
        for line in p.stdout:
            print(line.decode('utf-8').rstrip())

class Alignment(Atoms):
    def __init__(
        self,
        atoms,
        weights = None,
        biasing = False,
        bias_tol = 0.2,
        bias_scale = 1.e3,
        iteration = False,
        testing = False,
        records = 1,
        count = 10,
        trials = None,
    ):
        if isinstance(atoms, Atoms):
            self.__dict__.update(atoms.__dict__)
        else:
            raise TypeError('An Atoms object was expected as argument')
        if type(biasing) is not bool:
            raise TypeError('"biasing" must be a boolean')
        if type(iteration) is not bool:
            raise TypeError('"iteration" must be a boolean')
        if type(testing) is not bool:
            raise TypeError('"testing" must be a boolean')
        if type(records) is not int:
            raise TypeError('"records" must be an integer')
        if type(count) is not int:
            raise TypeError('"count" must be an integer')
        if type(bias_tol) is not float:
            raise TypeError('"bias_tol" must be a real')
        if type(bias_scale) is not float:
            raise TypeError('"bias_scale" must be a real')
        if trials is None:
            molalignlib.set_abort_flag(False)
        elif type(trials) is int:
            molalignlib.set_abort_flag(True)
            molalignlib.set_max_trials(trials)
        else:
            raise TypeError('"trials" must be an integer')
        if weights is None:
            self.weights = np.ones(len(atoms), dtype=float),
        elif type(weights) is not np.ndarray:
            self.weights = weights
        else:
            raise TypeError('"weights" must be a numpy array')
        self.records = records
        molalignlib.set_bias_flag(biasing)
        molalignlib.set_bias_scale(bias_scale)
        molalignlib.set_bias_tol(bias_tol)
        molalignlib.set_conv_flag(iteration)
        molalignlib.set_test_flag(testing)
        molalignlib.set_max_count(count)
    def sorted(self, other):
        if not isinstance(other, Atoms):
            raise TypeError('An Atoms object was expected as argument')
        if len(other) != len(self):
            raise ValueError('Argument does no have the right length')
        znums0 = self.get_atomic_numbers()
        znums1 = other.get_atomic_numbers()
        types0 = np.ones(len(self), dtype=int)
        types1 = np.ones(len(other), dtype=int)
        coords0 = self.get_positions().transpose()
        coords1 = other.get_positions().transpose()
        n, maplist, mapcount, mindist = molalignlib.sort_atoms(znums0, znums1, \
            types0, types1, coords0, coords1, self.weights, self.records)
        return [i - 1 for i in maplist.transpose()[:n]], mapcount[:n], mindist[:n]
    def aligned(self, other):
        if not isinstance(other, Atoms):
            raise TypeError('An Atoms object was expected as first argument')
        if len(other) != len(self):
            raise ValueError('First argument does no have the right length')
        znums0 = self.get_atomic_numbers()
        znums1 = other.get_atomic_numbers()
        types0 = np.ones(len(self), dtype=int)
        types1 = np.ones(len(other), dtype=int)
        coords0 = self.get_positions().transpose()
        coords1 = other.get_positions().transpose()
        travec, rotmat = molalignlib.align_atoms(znums0, znums1, types0, types1, \
             coords0, coords1, self.weights)
        coords1 = np.matmul(rotmat, coords1).transpose() + travec
        return Atoms(numbers=znums1, positions=coords1)
