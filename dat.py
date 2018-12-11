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

# from ovito.io import import_file as ovito_import_file
# from ovito.modifiers import AtomicStrainModifier

class Dat:
	"""
	.dat files are designed to save trajectory (position and velocity) data for
	simulations in 2D, with constant number N of particles. These are binary
	files, containing numbers organised according to the following scheme:
	|         FRAME (= CALL TO active_particles.dat.Dat.dump) 0         | ...
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
		variable : string
			Name of variable.

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
		(see active_particles.dat.Dat.variable)

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
		(see active_particles.dat.Dat.variable)

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

	def displacement(self, time0, time1, *particle):
		"""
		Returns array of displacement between times 'time0' and 'time1'.
		(see active_particles.dat.Dat.position)

		Parameters
		----------
		time0 : int
			Initial frame index.
		time1 : int
			Final frame index.

		Optional positional arguments
		-----------------------------
		particle : int
			Particle index.
			When called with particle indexes, function returns array of
			particles' displacement between frames 'time0' and 'time1' in the
			same order.

		Returns
		-------
		arr : self.element_type packing format Numpy array
			Array of displacement between frames 'time0' and 'time1'.
		"""

		return self.position(time1, *particle)\
			- self.position(time0, *particle)

class Gsd(HOOMDTrajectory):
	"""
	This class adds methods to the gsd.hoomd.HOOMDTrajectory class which reads
	.gsd trajectory file.
	"""

	def __init__(self, file, prep_frames=0, dimensions=2):
		"""
		Parameters
		----------
		file : file object
			Trajectory file. (.gsd)
		prep_frames : int
			Number of frames to ignore at beginning of .gsd file. (default: 0)
		dimensions : int
			Dimension of space. (default: 2)
		"""

		self.filename = file.name	# file name
		self.file = GSDFile(file)	# gsd file
		super().__init__(self.file)	# initialising gsd.hoomd.HOOMDTrajectory

		self.prep_frames = prep_frames
		self.dimensions = dimensions

		self.node = ovito_import_file(self.filename)	# OVITO ObjectNode used for nonaffine squared displacement computation

	def __getitem__(self, key):
		"""
		Parameters
		----------
		key : int or slice
			Index of trajectory frame.
			NOTE: To time is added the number of preparation frames
			self.prep_frames.

		Returns
		-------
		snapshot : gsd.hoomd.Snapshot object
			Snapshot(s) at time(s) key.
		"""

		if isinstance(key, slice):
			return super().__getitem__(slice(
				int(key.start + self.prep_frames) if key.start!=None else None,
				int(key.stop + self.prep_frames) if key.stop!=None else None,
				key.step))
		return super().__getitem__(int(key + self.prep_frames))

	def position(self, time, *particle, **kwargs):
		"""
		Returns array of position at frame 'time'.

		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particle : int
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

		positions = self[time].particles.position[:, :self.dimensions]	# positions at frame time
		if particle != ():												# consider only particles in particles
			positions = np.array(itemgetter(*particle)(positions))

		if 'centre' in kwargs:
			return relative_positions(positions, kwargs['centre'],
				self.box_size(time))	# positions with centre as centre
		return positions

	def velocity(self, time, *particle):
		"""
		Returns array of velocity at frame 'time'.

		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particle : int
			Particles indexes.
			When called with particles indexes, function returns array of
			particles' position at frame 'time' in the same order.

		Returns
		-------
		velocities : float Numpy array
			Array of velocities at frame 'time'.
		"""

		velocities = self[time].particles.velocity[:, :self.dimensions]		# velocities at frame time
		if particle == ():	return velocities								# returns all velocities
		return np.array(itemgetter(*particle)(velocities))					# velocities at frame time

	def diameter(self, time, *particle):
		"""
		Returns array of diameter at frame 'time'.

		Parameters
		----------
		time : int
			Frame index.

		Optional positional arguments
		-----------------------------
		particle : int
			Particles indexes.
			When called with particles indexes, function returns array of
			particles' diameters at frame 'time' in the same order.

		Returns
		-------
		diameters : float Numpy array
			Array of diameters at frame 'time'.
		"""

		diameters = self[time].particles.diameter			# diameters at frame time
		if particle == ():	return diameters				# returns all diameters
		return np.array(itemgetter(*particle)(diameters))	# diameters at frame time

	def box_size(self, time=0):
		"""
		Returns length of system box in first direction at time 'time'.

		Parameters
		----------
		time : int
			Frame index. (default: 0)

		Returns
		-------
		L : float
			Length of system box in fist direction.
		"""

		return self[time].configuration.box[0]

	def d2min(self, time0, time1, *particle):
		"""
		Returns nonaffine squared displacement computed by OVITO (see
		https://ovito.org/manual/particles.modifiers.atomic_strain.html and
		https://ovito.org/manual/python/modules/ovito_modifiers.html) between
		frames 'time0' and 'time1'.

		Parameters
		----------
		time0 : int
			Initial frame index.
		time1 : int
			Final frame index.

		Optional positional arguments
		-----------------------------
		particle : int
			Particles indexes.
			When called with particles indexes, function returns array of
			particles' diameters at frame 'time' in the same order.

		Returns
		-------
		d2min : float Numpy array
			Array of nonaffine square displacements between frames 'time0' and
			'time1'.
		"""

		self.node.modifiers.clear()											# clear modification pipeline
		self.node.modifiers.append(
			AtomicStrainModifier(
				output_nonaffine_squared_displacements=True,
				reference_frame=self.prep_frames + time0))					# add AtomicStrainModifier modifier to modification pipeline
		self.node.modifiers[-1].reference.load(self.filename)				# load trajectory file as reference
		self.node_out = self.node.compute(frame=self.prep_frames + time1)	# compute d2min

		d2min = self.node_out['Nonaffine Squared Displacement'].array	# array of nonaffine squared displacement
		if particle == ():	return d2min								# returns all nonaffine squared displacement
		return np.array(itemgetter(*particle)(d2min))					# non affine square displacements
