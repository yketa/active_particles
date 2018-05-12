"""
Module dat defines the object Dat, which allows one to write to and read from
.dat files.
"""

import numpy as np
import struct

class Dat:
	"""
	.dat files are designed to save trajectory (position and velocity) data for
	simulations in 2D, with constant number N of particles. These are binary
	files, containing numbers organised according to the following scheme:
	|                    FRAME (= CALL TO Dat.dump) 0                   | ...
	|           POSITIONS             |           VELOCITIES            | ...
	| PARTICLE 0 | ... | PARTICLE N-1 | PARTICLE 0 | ... | PARTICLE N-1 | ...
	|   x  |  y  | ... |   x   |  y   |   x  |  y  | ... |   x   |  y   | ...
	"""

	def __init__(self, data_file, N, element_type='d'):
		"""
		Parameters
		----------
		data_file : file object
			.dat file
			NOTE: In order to be able to both read and write, it is advised to
			open data_file in 'r+b' mode.
		N : int
			Number of particles.
		element_type : packing format
			Data packing format. (default: double float)
		"""

		self.file = data_file									# .dat file
		self.N = N												# number of particles
		self.element_type = element_type						# data picking format
		self.bytes_per_element = struct.calcsize(element_type)	# element_type number of bytes

	def dump(self, positions, velocities):
		"""
		Dump to file following the .dat file format (trajectory file).

		NOTE: self.file has to be open in 'wb' or 'r+b' mode.

		Parameters
		----------
		positions : (self.N, 2) shaped array like
			List of position coordinates.
		velocities : (self.N, 2) shaped array like
			List of velocity coordinates.
		"""

		for data in [positions, velocities]:							# for each variable to dump
			for particle in range(self.N):								# for each particle
				for coord in range(2):									# for each dimension of space
						self.file.write(struct.pack(self.element_type,	# dump to file
							data[particle][coord]))

	def get_value(self, time, particle, axis, variable='position'):
		"""
		Returns the projection on axis 'axis' of the variable 'variable' of
		particle 'particle' at the frame 'time'.

		NOTE: self.file has to be open in 'rb' or 'r+b' mode.

		Parameters
		----------
		time : int
			Frame index.
		particle : int
			Particle index.
		axis : int (either 0 for x-axis or 1 for y-axis)
			Axis index.
		variable : string (either 'position' or 'velocity')
			Name of the variable. (default: position)

		Returns
		-------
		val : self.element_type packing format
			Variable.
		"""

		try:
			inc_var = {'position':0, 'velocity':2*N}[variable]	# increment in bytes_per_element to accesss variable
		except KeyError:										# if variable is not 'position' or 'velocity'
			inc_var = 0											# consider variable as 'position'

		self.file.seek(self.bytes_per_element*(
			4*self.N*time + inc_var + 2*particle + axis))	# set file's current position according to frame, variable, number of particles, and axis
		return struct.unpack(self.element_type,
			self.file.read(self.bytes_per_element))[0]		# variable

	def get_array(self, time, variable='position'):
		"""
		Returns the (self.N, 2) array of the variable 'variable' for each
		particle at frame 'time'.

		NOTE: self.file has to be open in 'rb' or 'r+b' mode.

		Parameters
		----------
		time : int
			Frame index.
		variable : string (either 'position' or 'velocity')
			Name of the variable. (default: position)

		Returns
		-------
		arr : self.element_type packing format (self.N, 2) Numpy array
			Array of variable at frame 'time'.
		"""

		return np.reshape(list(map(
			lambda particle: list(map(
			lambda axis: self.get_value(time, particle, axis,
			variable=variable),
			range(2))),
			range(self.N))),
			(self.N, 2))

	def variable(self, time, *particle, variable='position'):
		"""
		Returns array of variable at frame 'time'.

		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particle : int
			Particle index.
			When called with particle indexes, function returns array of
			particles' variable at frame 'time' in the same order.

		Returns
		-------
		arr : self.element_type packing format Numpy array
			Array of variable at frame 'time'.
		"""

		if particle == (): return self.get_array(time, variable=variable)	# no particular indexes requested

		return np.reshape(list(map(
			lambda particle: list(map(
			lambda axis: self.get_value(time, particle, axis,
			variable=variable),
			range(2))),
			particle)),
			(len(particle), 2))	# variable at frame 'time' for particles 'particle'

	def position(self, time, *particle):
		"""
		Returns array of position at frame 'time'.
		(see dat.Dat.variable)

		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particle : int
			Particle index.
			When called with particle indexes, function returns array of
			particles' position at frame 'time' in the same order.

		Returns
		-------
		arr : self.element_type packing format Numpy array
			Array of position at frame 'time'.
		"""

		self.variable(time, *particle, variable='position')

	def velocity(self, time, *particle):
		"""
		Returns array of velocity at frame 'time'.
		(see dat.Dat.variable)

		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particle : int
			Particle index.
			When called with particle indexes, function returns array of
			particles' velocity at frame 'time' in the same order.

		Returns
		-------
		arr : self.element_type packing format Numpy array
			Array of velocity at frame 'time'.
		"""

		self.variable(time, *particle, variable='velocity')
