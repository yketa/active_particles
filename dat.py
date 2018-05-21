"""
Module dat defines the object Dat, which allows one to write to and read from
.dat files.
"""

import numpy as np
import struct
from operator import itemgetter

from active_particles.maths import relative_positions

from gsd.pygsd import GSDFile
from gsd.hoomd import HOOMDTrajectory

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
		self.N = int(N)											# number of particles
		self.element_type = element_type						# data picking format
		self.bytes_per_element = struct.calcsize(element_type)	# element_type number of bytes
		self.inc_var = {'position':0, 'velocity':2*self.N}		# increment in bytes_per_element to accesss variable

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

	def get_value(self, time, particle, axis, inc_var):
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
		inc_var : int
			Increment in bytes_per_element to accesss variable.

		Returns
		-------
		val : self.element_type packing format
			Variable.
		"""

		self.file.seek(self.bytes_per_element*(
			4*self.N*time + inc_var + 2*particle + axis))	# set file's current position according to frame, variable, number of particles, and axis
		return struct.unpack(self.element_type,
			self.file.read(self.bytes_per_element))[0]		# variable

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

		if particle == (): particle = range(self.N)	# no particular indexes requested

		inc_var = self.inc_var[variable]	# increment in bytes_per_element to accesss variable

		return np.reshape(list(map(
			lambda particle: list(map(
			lambda axis: self.get_value(time, particle, axis, inc_var),
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

		return self.variable(time, *particle, variable='position')

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

		return self.variable(time, *particle, variable='velocity')

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

		self.file = GSDFile(file)	# gsd file
		super().__init__(self.file)	# initialising gsd.hoomd.HOOMDTrajectory
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

		positions = self[time].particles.position[:, :self.dimensions]	# positions at frame time
		if particles != ():												# consider only particles in particles
			positions = np.array(itemgetter(*particles)(positions))

		if 'centre' in kwargs:
			box_dim = self[time].configuration.box[0]						# box dimensions
			return relative_positions(positions, kwargs['centre'], box_dim)	# positions with centre as centre
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

		velocities = self[time].particles.velocity[:, :self.dimensions])	# velocities at frame time
		if particles == ():	return velocities								# returns all positions
		return np.array(itemgetter(*particles)(velocities))					# velocities at frame time
