#!/usr/bin/env python

import math
import numpy as np


def _angle_between(a, b):
	a = np.copy(np.array(a))
	b = np.copy(np.array(b))
	if np.linalg.norm(a) == 0 or np.linalg.norm(b) == 0:
		raise ValueError('Got zero vector.')
	a /= np.linalg.norm(a)
	b /= np.linalg.norm(b)
	angle = np.arccos(np.dot(a, b))

	if np.isnan(angle):
		if np.isnan(np.min(a)) or np.isnan(np.min(b)):
			raise ValueError('NaN in vectors.')
		if (a == b).all():
			return 0.0
		else:
			return np.pi
	return angle


def hmatrix_to_abc(h_matrix, degrees=False):
	result = np.zeros(6)
	for i in range(3):
		result[i] = np.linalg.norm(h_matrix[:, i])
	result[3] = _angle_between(h_matrix[:, 1], h_matrix[:, 2])
	result[4] = _angle_between(h_matrix[:, 0], h_matrix[:, 2])
	result[5] = _angle_between(h_matrix[:, 0], h_matrix[:, 1])
	if degrees:
		result[3:] = map(math.degrees, result[3:])
	return result


def abc_to_hmatrix(a, b, c, alpha, beta, gamma, degrees=True):
	if degrees:
		alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
	result = np.zeros((3, 3))

	a = np.array((a, 0, 0))
	b = b * np.array((math.cos(gamma), math.sin(gamma), 0))
	bracket = (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
	c = c * np.array((math.cos(beta), bracket, math.sqrt(math.sin(beta) ** 2 - bracket ** 2)))

	result[:, 0] = a
	result[:, 1] = b
	result[:, 2] = c

	return result


def repeat_vector(h_matrix, repeat_a, repeat_b, repeat_c):
	return h_matrix[:, 0] * repeat_a + h_matrix[:, 1] * repeat_b + h_matrix[:, 2] * repeat_c


def box_vertices(h_matrix, repeat_a, repeat_b, repeat_c):
	vertices = np.zeros((8, 3))
	vertices[0, :] = repeat_vector(h_matrix, repeat_a, repeat_b, repeat_c)
	vertices[1, :] = repeat_vector(h_matrix, repeat_a, repeat_b, repeat_c + 1)
	vertices[2, :] = repeat_vector(h_matrix, repeat_a, repeat_b + 1, repeat_c)
	vertices[3, :] = repeat_vector(h_matrix, repeat_a, repeat_b + 1, repeat_c + 1)
	vertices[4, :] = repeat_vector(h_matrix, repeat_a + 1, repeat_b + 1, repeat_c + 1)
	vertices[5, :] = repeat_vector(h_matrix, repeat_a + 1, repeat_b + 1, repeat_c)
	vertices[6, :] = repeat_vector(h_matrix, repeat_a + 1, repeat_b, repeat_c + 1)
	vertices[7, :] = repeat_vector(h_matrix, repeat_a + 1, repeat_b, repeat_c)
	return vertices


def cell_volume(h_matrix):
	ab = np.cross(h_matrix[:, 0], h_matrix[:, 1])
	return np.abs(np.dot(ab, h_matrix[:, 2]))


def cartesian_to_scaled_coordinates(coordinates, h_matrix):
	h = np.linalg.inv(h_matrix)
	for i in range(len(coordinates)):
		coordinates[i] = (h * coordinates[i]).sum(axis=1)
	return coordinates


def scaled_to_cartesian_coordinates(coordinates, h_matrix):
	for i in range(len(coordinates)):
		coordinates[i] = (h_matrix * coordinates[i]).sum(axis=1)
	return coordinates


def cell_multiply(coord, x, y, z, h_matrix=None, scaling_in=False, scaling_out=False):
	for i in (x, y, z):
		if i < 1 or int(i) != i:
			raise ValueError('Invalid image count.')

	# copy input
	coord = np.copy(coord)

	# prepare data
	factor = x * y * z
	atoms = coord.shape[0]
	newcoord = np.zeros((atoms * factor, 3))
	offset = 1

	# copy data
	if scaling_out:
		if not scaling_in:
			if h_matrix is None:
				raise TypeError('H matrix has to be given for partially cartesian data.')
			coord = cartesian_to_scaled_coordinates(coord, h_matrix)
		newcoord[:atoms] = coord

		for i_x in range(x):
			for i_y in range(y):
				for i_z in range(z):
					if i_x == i_y == i_z == 0:
						continue
					newcoord[atoms * offset:atoms * (offset + 1)] = coord
					newcoord[atoms * offset:atoms * (offset + 1), 0] += (i_x)
					newcoord[atoms * offset:atoms * (offset + 1), 1] += (i_y)
					newcoord[atoms * offset:atoms * (offset + 1), 2] += (i_z)
					offset += 1
		newcoord[:, 0] /= x
		newcoord[:, 1] /= y
		newcoord[:, 2] /= z
	else:
		if h_matrix is None:
			raise TypeError('H matrix has to be given for cartesian data.')
		if scaling_in:
			coord = scaled_to_cartesian_coordinates(coord, h_matrix)
		newcoord[:atoms] = coord
		for i_x in range(x):
			for i_y in range(y):
				for i_z in range(z):
					if i_x == i_y == i_z == 0:
						continue
					vector = h_matrix[:, 0] * i_x + h_matrix[:, 1] * i_y + h_matrix[:, 2] * i_z
					newcoord[atoms * offset:atoms * (offset + 1)] = coord + vector
					offset += 1

	return newcoord