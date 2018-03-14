#! /home/yketa/miniconda3/bin/python3.6

import struct

def getval(data, N, time, particle, axis, variable='position', element_type='d', bytes_per_element=8):
	"""
	Returns the projection on axis 'axis' of the variable 'variable' of particle 'particle' at the frame 'time' for an opened binary file 'data' corresponding to a simulation with 'N' particles.
	By default, each value is a float saved as a double (element_type='d') which then takes bytes_per_element=8 bytes.
	"""

	data.seek(bytes_per_element*(4*N*time + {'position':0, 'velocity':2*N}[variable] + 2*particle + axis)) # set the file's position
	return struct.unpack(element_type, data.read(bytes_per_element))[0]
