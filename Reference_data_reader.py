import scipy.io
import numpy as np
mat = scipy.io.loadmat('Reference_data.mat') # this creates a dictionary with the data

#print(mat.keys())
Flight_Data = mat["flightdata"] #numpy.ndarray
print(Flight_Data[0])
