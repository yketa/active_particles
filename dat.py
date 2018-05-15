"""
Module dat defines the object Dat, which allows one to write to and read from
.dat files.
"""

import numpy as np
import struct
from collections import OrderedDict
from operator import itemgetter

from gsd.hoomd import HOOMDTrajectory

class Dat:
	"""
	.dat files are designed to save time-dependent variables (e.g., position
	and velocity) for systems of any dimensions, with constant number N of
	particles. These are binary files, containing numbers organised according
	to the following scheme:

	|          FRAME (=CALL TO active_particles.dat.Dat.dump) 1          |..
	|           VARIABLE 1           |..|           VARIABLE M           |..
	|  PARTICLE 1  |..|  PARTICLE N  |..|  PARTICLE 1  |..|  PARTICLE N  |..
	|x_1|..|x_{d_1}|..|x_1|..|x_{d_1}|..|x_1|..|x_{d_M}|..|x_1|..|x_{d_M}|..

	with M the number of variables, N the number of particles, d_i the
	dimension of the i_th variable.
	"""

	def __init__(self, data_file, N, variables=OrderedDict([('position', 2),
		('velocity', 2)]), element_type='d'):
		"""
		Parameters
		----------
		data_file : file object
			.dat file
			NOTE: In order to be able to both read and write, it is advised to
			open data_file in 'r+b' mode.
		N : int
			Number of particles.
		variables : ordered dictionary
			Variable names as keys and dimensions as values. (default:
			{'position': 2, 'velocity': 2})
			NOTE: It is strongly advised to have 'position' as first variable.
		element_type : packing format
			Data packing format. (default: double float)
		"""

		self.file = data_file									# .dat file
		self.N = int(N)											# number of particles
		self.variables = variables								# variables dictionary
		self.v_number = len(self.variables)						# number of variables
		self.v_names = list(self.variables.keys())				# variable names
		self.v_dim = list(self.variables.values())				# variable dimensions
		self.element_type = element_type						# data picking format
		self.bytes_per_element = struct.calcsize(element_type)	# element_type number of bytes

	def dump(self, *var_data):
		"""
		Dump a frame to file following the .dat file format.

		NOTE: self.file has to be open in 'wb' or 'r+b' mode.

		Positional arguments
		--------------------
		var_data : (self.N, self.variables[variable name]) shaped array like
			Array of variable coordinates.
			NOTE: Arrays of variables have to be passed in the order of
			self.variables.
		"""

		for var in range(self.v_number):										# for all variables to dump
			data = np.reshape(var_data[var], (self.N, self.v_dim[var]))			# reshaping variable according to known dimension
			for particle in range(self.N):										# for all particles
				for coord in range(self.v_dim[var]):							# for each dimension of variable
					self.file.write(struct.pack(self.element_type,				# dump to file
						data[particle, coord]))

	def get_value(self, time, particle, axis=0, variable='position'):
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
		axis : int
			Axis index. (default: 0)
		variable : string (either 'position' or 'velocity')
			Name of the variable. (default: 'position')

		Returns
		-------
		val : self.element_type packing format
			Variable.
		"""

		try:
			index_var = self.v_names.index(variable)	# index of variable in variables ordered dictionary
		except KeyError:								# if variable is not not known
			index_var = 0								# consider variable as first variable

		axis = 0 if axis > self.v_dim[index_var] else axis	# set axis as 0 if input axis is greater than variable dimension

		self.file.seek(self.bytes_per_element*(					# set file's current position
			time*np.sum(self.v_dim, dtype=int)*self.N			# to desired frame
			+ np.sum(self.v_dim[:index_var], dtype=int)*self.N	# to desired variable
			+ self.v_dim[index_var]*particle					# to desired particle
			+ axis))											# to desired axis

		return struct.unpack(self.element_type,
			self.file.read(self.bytes_per_element))[0]		# variable

	def get_value_vec(self, time, particles=None, axes=None,
		variable='position'):
		"""
		Returns the projection on axes in 'axis' of the variable 'variable' of
		particles in 'particle' at the frame 'time'.

		NOTE: self.file has to be open in 'rb' or 'r+b' mode.

		Parameters
		----------
		time : int
			Frame index.
		particles : int array-like
			Particles index array. (default: None)
			None or empty tuple is equivalent to all particles (particles =
			range(self.N)).
			DEFAULT: None
		axes : int array-like
			Axes index array. (default: None)
			None or empty tuple is equivalent to all axes (axes =
			range(self.variables[variable]))
		variable : string
			Name of the variable. (default: 'position')

		Returns
		-------
		val : self.element_type packing format (len(particles), len(axes))
		Numpy array
			List of variable.
		"""

		if particles == None or particles == (): particles = range(self.N)
		if axes == None or axes == ():
			try:
				axes = range(self.variables[variable])
			except:			# if variable is not not known
				axes = (0)	# return first axis value


		return np.reshape(list(map(
			lambda particle: list(map(
			lambda axis: self.get_value(time, particle, axis=axis,
			variable=variable),
			axes)),
			particles)),
			(len(particles), len(axes)))

	def get_array(self, time, variable='position'):
		"""
		Returns the (self.N, self.variables[variable]) array of the variable
		'variable' for each particle at frame 'time'.
		(see active_particles.dat.Dat.get_value_vec)

		NOTE: self.file has to be open in 'rb' or 'r+b' mode.

		Parameters
		----------
		time : int
			Frame index.
		variable : string
			Name of the variable. (default: 'position')

		Returns
		-------
		arr : self.element_type packing format
		(self.N, self.variables[variable]) Numpy array
			Array of variable at frame 'time'.
		"""

		return self.get_value_vec(time, variable=variable)

	def variable(self, time, *particles, variable='position'):
		"""
		Returns array of variable at frame 'time'.
		(see active_particles.dat.Dat.get_value_vec)


		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particles : int
			Particles indexes.
			When called with particles indexes, function returns array of
			particles' variable at frame 'time' in the same order.

		Returns
		-------
		arr : self.element_type packing format Numpy array
			Array of variable at frame 'time'.
		"""

		return self.get_value_vec(time, particles=particles, variable=variable)

	def position(self, time, *particles):
		"""
		Returns array of position at frame 'time'.
		(see active_particles.dat.Dat.variable)

		NOTE: 'position' has to be in self.variables

		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particles : int
			Particles indexes.
			When called with particles indexes, function returns array of
			particles' position at frame 'time' in the same order.

		Returns
		-------
		arr : self.element_type packing format Numpy array
			Array of position at frame 'time'.
		"""

		return self.variable(time, *particles, variable='position')

	def velocity(self, time, *particles):
		"""
		Returns array of velocity at frame 'time'.
		(see dat.Dat.variable)

		NOTE: 'velocity' has to be in self.variables

		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particles : int
			Particles indexes.
			When called with particles indexes, function returns array of
			particles' velocity at frame 'time' in the same order.

		Returns
		-------
		arr : self.element_type packing format Numpy array
			Array of velocity at frame 'time'.
		"""

		return self.variable(time, *particles, variable='velocity')

class Gsd(HOOMDTrajectory):
	"""
	This class adds methods to the gsd.hoomd.HOOMDTrajectory class which reads
	.gsd trajectory file.
	"""

	def __init__(self, file, dimensions=2):
		"""
		Parameters
		----------
		file : gsd.pygsd.GSDFile object
			Trajectory file. (.gsd)
		dimensions : int
			Dimension of space. (default: 2)
		"""

		super().__init__(file)	# initialising gsd.hoomd.HOOMDTrajectory
		self.dimensions = dimensions

	def position(self, time, *particles, **kwargs):
		"""
		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particles : int
			Particles indexes.
			When called with particles indexes, function returns array of
			particles' position at frame 'time' in the same order.

		Optional keyword arguments
		--------------------------
		centre : (self.dimensions,) array
			Define new centre position.

		Returns
		-------
		positions : float Numpy array
			Array of positions at frame 'time'.
		"""

		time = int(time)	# avoids crash when calling self.__getitem__

		if particles == ():	particles = range(self[time].particles.N)	# returns all positions
		positions = itemgetter(*particles)(
			self[time].particles.position[:, :self.dimensions])			# positions at frame time

		if 'centre' in kwargs:
			box_dim = self[time].configuration.box[:self.densions]	# box dimensions
			return (positions - np.array(kwargs['centre'])
				+ box_dim/2)%box_dim - box_dim/2					# positions with centre as centre
		return positions

	def velocity(self, time, *particles):
		"""
		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particles : int
			Particles indexes.
			When called with particles indexes, function returns array of
			particles' position at frame 'time' in the same order.

		Returns
		-------
		velocities : float Numpy array
			Array of velocities at frame 'time'.
		"""

		time = int(time)	# avoids crash when calling self.__getitem__

		if particles == ():	particles = range(self[time].particles.N)	# returns all positions
		return itemgetter(*particles)(
			self[time].particles.velocity[:, :self.dimensions])			# velocities at frame time
