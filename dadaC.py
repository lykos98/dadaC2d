import ctypes as ct
import numpy as np  
import os

ctFloatType = ct.c_double
ctIdxType = ct.c_uint64

class HeapNode(ct.Structure):
    _fields_ = [
        ("value", ctFloatType),
        ("array_idx", ctIdxType)
    ]


class Heap(ct.Structure):
    _fields_ = [
        ("N", ctIdxType),
        ("count", ctIdxType),
        ("data", ct.POINTER(HeapNode))
    ]

class luDynamicArray(ct.Structure):
    _fields_ = [
        ("data", ct.POINTER(ctIdxType)),
        ("size", ctIdxType),
        ("count", ctIdxType)
    ]

class DatapointInfo(ct.Structure):
    _fields_ = [
        ("g", ctFloatType),
        ("ngbh", Heap),
        ("array_idx", ctIdxType),
        ("log_rho", ctFloatType),
        ("log_rho_c", ctFloatType),
        ("log_rho_err", ctFloatType),
        ("kstar", ctIdxType),
        ("is_center", ct.c_int),
        ("cluster_idx", ct.c_int)
    ]

    def __repr__(self):
        return f"{{ g: {self.g}, ngbh: {self.ngbh}, array_idx: {self.array_idx}, log_rho: {self.log_rho}, log_rho_c: {self.log_rho_c}, log_rho_err: {self.log_rho_err}, kstar: {self.kstar}, isCenter: {self.is_center}, clusterIdx: {self.cluster_idx} }}"

class Border_t(ct.Structure):
    _fields_ = [
        ("idx", ctIdxType),
        ("density", ctFloatType),
        ("error", ctFloatType)
    ]

class SparseBorder_t(ct.Structure):
    _fields_ = [
        ("i", ctIdxType),
        ("j", ctIdxType),
        ("idx", ctIdxType),
        ("density", ctFloatType),
        ("error", ctFloatType)
    ]

class AdjList(ct.Structure):
    _fields_ = [
        ("count", ctIdxType),
        ("size", ctIdxType),
        ("data", ct.POINTER(SparseBorder_t))
    ]

class Clusters(ct.Structure):
    _fields_ = [
        ("UseSparseBorders", ct.c_int),
        ("SparseBorders", ct.POINTER(AdjList)),
        ("centers", luDynamicArray),
        ("borders", ct.POINTER(ct.POINTER(Border_t))),
        ("__borders_data", ct.POINTER(Border_t)),
        ("n", ctIdxType)
    ]

class Data():
    def __init__(self, data : np.array):
        """Data object for dadaC library

        Args:
            data (np.array): 2d array to use in searching for clusters 

        Raises:
            TypeError: Raises TypeError if a type different from a matrix is passed 
        """
        path = os.path.join(os.path.dirname(__file__), "bin/libdadac.so")
        self.lib = ct.CDLL(path)

        global ctFloatType, ctIdxType
        s = self.lib.FloatAndUintSize()

        self.__useFloat32 = s < 1
        self.__useInt32   = (s == 0) or (s == 2)

        #initialize class
        self.state = {
                    "ngbh" : False,
                    "id" : False,
                    "density" : False,
                    "clustering" : False,
                    "useFloat32": self.__useFloat32,
                    "useInt32": self.__useInt32,
                    "useSparse": None,
                    "computeHalo": None 
                    }
        if self.__useFloat32:
            self.data = data.astype(np.float32)
            ctFloatType = ct.c_float
        else:
            self.data = data.astype(np.float64)

        if self.__useInt32:
            ctIdxType = ct.c_uint32

        #retrieve function pointers form .so file

        self.__NgbhSearch = self.lib.NgbhSearch
        self.__NgbhSearch.argtypes = [np.ctypeslib.ndpointer(ctFloatType), ct.c_uint64, ct.c_uint64, ct.c_uint64 ]
        self.__NgbhSearch.restype  = ct.POINTER(DatapointInfo)

        self.__idEstimate = self.lib.idEstimate
        self.__idEstimate.argtypes = [ct.POINTER(DatapointInfo), ct.c_uint64]
        self.__idEstimate.restype  =  ct.c_double


        self.__computeRho = self.lib.computeRho
        self.__computeRho.argtypes = [ct.POINTER(DatapointInfo), ct.c_double, ct.c_uint64]

        self.__computeCorrection = self.lib.computeCorrection
        self.__computeCorrection.argtypes = [ct.POINTER(DatapointInfo), ctIdxType, ct.c_double]

        self.__H1 = self.lib.Heuristic1
        self.__H1.argtypes = [ct.POINTER(DatapointInfo), np.ctypeslib.ndpointer(ctFloatType), ct.c_uint64]
        self.__H1.restype = Clusters

        self.__ClustersAllocate = self.lib.Clusters_allocate
        self.__ClustersAllocate.argtypes = [ct.POINTER(Clusters), ct.c_int]

        self.__H2 = self.lib.Heuristic2
        self.__H2.argtypes = [ct.POINTER(Clusters), ct.POINTER(DatapointInfo)]

        self.__H3 = self.lib.Heuristic3
        self.__H3.argtypes = [ct.POINTER(Clusters), ct.POINTER(DatapointInfo), ct.c_double, ct.c_int]
        
        self.__freeDatapoints = self.lib.freeDatapointArray
        self.__freeDatapoints.argtypes = [ct.POINTER(DatapointInfo), ct.c_uint64]

        self.__freeClusters = self.lib.Clusters_free
        self.__freeClusters.argtypes = [ct.POINTER(Clusters)]


        if len(self.data.shape) != 2:
            raise TypeError("Please provide a 2d numpy array")


        self.__datapoints     = None
        self.__clusters       = None
        self.n              = self.data.shape[0]
        self.dims           = self.data.shape[1]
        self.k              = None
        self.id             = None

        self.clusterAssignment = None
        self.neighbors         = None
        self.borders           = None
        self.density           = None
        self.densityError      = None



    def computeNeighbors(self, k : int):
        
        """Compute the k nearest neighbors of each point

        Args:
            k (int): Number of neighbors to compute for each point 
        """
        self.k = k
        self.__datapoints = self.__NgbhSearch(self.data, self.n, self.dims, self.k)
        self.state["ngbh"] = True
        self.neighbors = None

    def computeIDtwoNN(self):

        """ Compute the intrinsic dimension of the dataset via the TWO Nearest Neighbors method.
            Ref. paper 

        Raises:
            ValueError: Raises value error if neighbors are not computed yet, use `Data.computeNeighbors()` 
        """

        if not self.state["ngbh"]:
            raise ValueError("Please compute Neighbors before calling this function")
        self.id = self.__idEstimate(self.__datapoints,self.n)
        self.state["id"] = True

    def computeDensity(self):

        """Compute density value for each point

        Raises:
            ValueError: Raises value error if ID is not computed, use `Data.computeIDtwoNN` method
        """
        if not self.state["id"]:
            raise ValueError("Please compute ID before calling this function")
        self.__computeRho(self.__datapoints, self.id, self.n)
        self.state["density"] = True
        self.density = None
        self.densityError = None

    def computeClusteringADP(self,Z : float, halo = True, useSparse = "auto"):

        """Compute clustering via the Advanced Density Peak method

        Args:
            Z (float): Z value for the method 
            useSparse (str): optional [``auto``,`yes`,`no`], use sparse implementation of border storage between clusters. Memory usage for big datsets is significant. 

        Raises:
            ValueError: Raises value error if density is not computed, use `Data.computeDensity()` method
        """
        if not self.state["density"]:
            raise ValueError("Please compute density before calling this function")
        if useSparse == "auto":
            if self.n > 2e6:
                self.state["useSparse"] = True
            else: 
                self.state["useSparse"] = False 
        elif useSparse == "y":
            self.state["useSparse"] = True
        else:
            self.state["useSparse"] = False

        self.state["computeHalo"] = halo 
        self.Z = Z
        self.__computeCorrection(self.__datapoints, self.n, self.Z)
        self.__clusters = self.__H1(self.__datapoints,self.data, self.n)
        self.__ClustersAllocate(ct.pointer(self.__clusters), 1 if self.state["useSparse"] else 0)
        self.__H2(ct.pointer(self.__clusters), self.__datapoints)
        self.__H3(ct.pointer(self.__clusters), self.__datapoints, self.Z, 1 if halo else 0 )
        self.state["clustering"] = True
        self.clusterAssignment = None

    def getClusterAssignment(self) -> list:


        """Retrieve cluster assignment

        Raises:
            ValueError: Raises error if clustering is not computed, use `Data.computeClusteringADP(Z)` 

        Returns:
            List of cluster labels
            
        """
        if self.clusterAssignment is None:
            if self.state["clustering"]:
                self.clusterAssignment = np.array([int(self.__datapoints[j].cluster_idx) for j in range(self.n)])
                return self.clusterAssignment
            else:
                raise ValueError("Clustering is not computed yet")
        else:
            return self.clusterAssignment

    def getBorders(self):
        raise NotImplemented("It's difficult I have to think about it")

    def getDensity(self) -> list:

        """Retrieve list of density values

        Raises:
            ValueError: Raise error if density is not computed, use `Data.computeDensity()` 

        Returns:
            List of density values
            
        """
        if self.density is None:
            if self.state["density"]:
                self.density = np.array([float(self.__datapoints[j].log_rho) for j in range(self.n)])
                return self.density
            else:
                raise ValueError("Density is not computed yet")
        else:
            return self.density

    def getDensityError(self):
        """Retrieve list of density error values

        Raises:
            ValueError: Raise error if density is not computed, use `Data.computeDensity()` 

        Returns:
            List of density error values
            
        """
        if self.densityError is None:
            if self.state["density"]:
                self.densityError = np.array([float(self.__datapoints[j].log_rho_err) for j in range(self.n)])
                return self.densityError
            else:
                raise ValueError("Density Error is not computed yet")
        else:
            return self.densityError
    
    def getNeighbors(self) -> list:
        """Retrieve k Nearest Neighbors of each point and their associated distance

        Raises:
            ValueError: Raise error if neighbors are not computed, use `Data.computeNeighbors(k)` 

        Returns:
            Returns lists of neighbors and distances            
        """
        if self.neighbors is None:
            if self.state["ngbh"]:
                self.neighbors = (
                            np.array([[int(self.__datapoints[j].ngbh.data[kk].array_idx) for kk in range(self.k)] for j in range(self.n)]),
                            np.array([[float(self.__datapoints[j].ngbh.data[kk].value)**(0.5) for kk in range(self.k)] for j in range(self.n)])
                        )
                return self.neighbors
            else:
                raise ValueError("Density is not computed yet")
        else:
            return self.neighbors

    def __del__(self):
        if not self.__datapoints is None:
            self.__freeDatapoints(self.__datapoints, self.n)
        if not self.__clusters is None:
            self.__freeClusters(ct.pointer(self.__clusters))
