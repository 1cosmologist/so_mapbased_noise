import os
try:
    OMP_NUM_THREADS=os.environ['OMP_NUM_THREADS']
except:
    print('OMP_NUM_THREADS not set')
import numba as nb
import numpy as np
from mpi4py import MPI
import sys

class MPIHandler:
    def __init__(self):
        """
        MPIHandler class initializer.

        Parameters
        ----------
        None

        Notes
        -----
        This class is the interface to the MPI world.

        Attributes
        ----------
        universe_comm : MPI communicator
            The MPI communicator for the entire universe of processes.
        universe_rank : int
            The rank of the local process in the universe of processes.
        universe_size : int
            The total number of processes in the universe of processes.
        """
        self.universe_comm = MPI.COMM_WORLD 
        self.universe_rank = self.universe_comm.Get_rank()
        self.universe_size = self.universe_comm.Get_size()
        self.nworlds       = self.universe_size

        if self.universe_rank == 0:
            print(f'Total processes: {self.universe_size}') ; sys.stdout.flush() 
            print(f'Number of threads: {nb.get_num_threads()}, {OMP_NUM_THREADS}'); sys.stdout.flush()

        return

    def split(self, procs_per_mpi, total_mcs):
        """
        Split the universe of processes into a specified number of worlds.

        Parameters
        ----------
        procs_per_mpi : int
            The number of processes to assign to each world.
        total_mcs : int
            The total number of worlds to create.

        Notes
        -----
        This function splits the universe of processes into a specified number of
        worlds. The number of worlds created is the minimum of the number of processes
        in the universe divided by the number of processes per world and the total
        number of worlds requested. The local process is then assigned to a world
        and the rank of the local process in its world is stored in the rank
        attribute. In addition, the total number of processes in the world is
        stored in the totProc attribute.

        The color attribute is the identity of the world
        if the local process is assigned to a world, otherwise it is set to MPI.UNDEFINED.
        The do_compute attribute is set to True if the local process is assigned
        to a world, otherwise it is set to False.
        """
        self.nworlds = int((self.universe_size-1) / procs_per_mpi)
        self.nworlds = min(self.nworlds, total_mcs)

        if self.universe_rank == 0:
            print(f'Total processes: {self.universe_size}') ; sys.stdout.flush() 
            print(f'Number of parallel MCs: {self.nworlds}') ; sys.stdout.flush()
            print(f'Proc per parallel MC: {procs_per_mpi}') ; sys.stdout.flush()
            print(f'Number of parallel computations: {self.nworlds}') ; sys.stdout.flush()

        self.color = MPI.UNDEFINED
        self.do_compute = (self.universe_rank > 0) and (self.universe_rank <= (self.nworlds * procs_per_mpi)) 
        if self.do_compute:
            self.color = int((self.universe_rank - 1) / procs_per_mpi)

        self.comm = self.universe_comm.Split(self.color)
        self.rank = MPI.UNDEFINED 
        self.totProc = MPI.UNDEFINED
        if self.do_compute:
            self.rank = self.comm.Get_rank()
            self.totProc = self.comm.Get_size()

            print(f'Rank={self.rank}; Number of threads: {nb.get_num_threads()}, {OMP_NUM_THREADS}'); sys.stdout.flush()

        return
    
    def divide_between_worlds(self, num2divide, starting_offset=0):
        """
        Divide a given number into slabs across the parallel worlds.

        Parameters
        ----------
        num2divide : int or np.int8 or np.int16 or np.int32 or np.int64
            The number to be divided into slabs.
        starting_offset : int, optional
            The starting offset for the slabs. Defaults to 0.

        Returns
        -------
        tuple
            A tuple of two arrays. The first array contains the starting index for each slab,
            and the second array contains the number of elements in each slab.

        Notes
        -----
        This function divides the given number into slabs across the parallel worlds.
        The number of slabs is determined by the number of parallel worlds, and the size of each slab
        is determined by the number of elements to be divided.
        """
        
        if self.nworlds == self.universe_size: print("WARNING: Worlds not split!")
        if not type(num2divide) in [int, np.int8, np.int16, np.int32, np.int64]:
            print("ERROR: Cannot divide a non-integer number")

        slab_min = num2divide // self.nworlds      # Q
        iter_remains = np.mod(num2divide, self.nworlds)    # R 

        #  SZ = R x (Q + 1) + (P-R) x Q 
        slab_per_Proc = np.zeros((self.nworlds,), dtype=np.int64)  # P = len(zslab_per_Proc)

        slab_per_Proc[0:int(iter_remains)] = slab_min + 1     # R procs together get (Q+1)xR z slabs
        slab_per_Proc[int(iter_remains):]  = slab_min 

        starts_per_Proc = np.zeros((self.nworlds,), dtype=np.int64)
        starts_per_Proc[0] = starting_offset
        starts_per_Proc[1:] = np.cumsum(slab_per_Proc)[:self.nworlds-1]

        return starts_per_Proc, slab_per_Proc
    
    def divide_between_procs(self, num2divide, starting_offset=0):
        """
        Divide a given number into slabs across the processes within a world.

        Parameters
        ----------
        num2divide : int or np.int8 or np.int16 or np.int32 or np.int64
            The number to be divided into slabs.
        starting_offset : int, optional
            The starting offset for the slabs. Defaults to 0.

        Returns
        -------
        tuple
            A tuple of two arrays. The first array contains the starting index for each slab,
            and the second array contains the number of elements in each slab.

        Notes
        -----
        This function divides the given number into slabs across the processes within a world.
        The number of slabs is determined by the number of processes within the world, and the size of each slab
        is determined by the number of elements to be divided.
        """
        if not isinstance(self.totProc, (int, np.int8, np.int16, np.int32, np.int64, MPI.INT64_T, MPI.INT32_T, MPI.INT16_T, MPI.INT8_T, MPI.UINT64_T, MPI.UINT32_T, MPI.UINT16_T, MPI.UINT8_T)):
            print("ERROR: Worlds not split! Aborting...")
            exit()

        if not type(num2divide) in [int, np.int8, np.int16, np.int32, np.int64]:
            print("ERROR: Cannot divide a non-integer number")

        slab_min = num2divide // self.totProc      # Q
        iter_remains = np.mod(num2divide, self.totProc)    # R 

        #  SZ = R x (Q + 1) + (P-R) x Q 
        slab_per_Proc = np.zeros((self.totProc,), dtype=np.int64)  # P = len(zslab_per_Proc)

        slab_per_Proc[0:int(iter_remains)] = slab_min + 1     # R procs together get (Q+1)xR z slabs
        slab_per_Proc[int(iter_remains):]  = slab_min 

        starts_per_Proc = np.zeros((self.totProc,), dtype=np.int64)
        starts_per_Proc[0] = starting_offset
        starts_per_Proc[1:] = np.cumsum(slab_per_Proc)[:self.totProc-1]

        return starts_per_Proc, slab_per_Proc