#!/usr/bin/env python

import numpy as np
import gzip
import argparse

parser = argparse.ArgumentParser(description='Reads CP2K PDOS files')
parser.add_argument('file', type=argparse.FileType('r'), help='PDOS output file to read.')

# range options
parser.add_argument('--limit', default=-1, type=int, help='The maximum number of frames to read from the PDOS files. Includes skipped frames. -1 = no limit.')
parser.add_argument('--skip', default=0, type=int, help='The number of frames to skip from the beginning of the PDOS files. 0 = none.')

# action options
parser.add_argument('--orbital_center', default=None, type=str, help='Calculate the occupied center for these orbitals.')
parser.add_argument('--orbital_cutoff', default=-1, type=float, help='Exclude orbitals with eigenvalues below this value. Unit: Hartree.')
parser.add_argument('--visualise', action='store_true', help='Whether to visualise the results with matplotlib.')
parser.add_argument('--relative_fermi', action='store_true', help='Whether to shift all values such that E_Fermi = 0.')
parser.add_argument('--relative_homo', action='store_true', help='Whether to shift all values such that E_HOMO = 0.')
parser.add_argument('--fermi_file', type=file, help='PDOS file to read the Fermi level from for each frame.')
parser.add_argument('--add', type=file, nargs='+', help='PDOS files to combine (same orbitals are merged)')
parser.add_argument('--subtract', type=file, nargs='+', help='PDOS files to combine with inverted sign (same orbitals are merged)')
parser.add_argument('--dos', action='store_true', help='Calculate the DOS averaged over a set of frames.')
parser.add_argument('--dos_cutoff', default=-1, type=float, help='Exclude orbitals with eigenvalues below this value. Unit: Hartree.')
parser.add_argument('--dos_smear', default=0.01, type=float, help='Gaussian smearing of the DOS. Unit: eV.')
parser.add_argument('--binwidth', default=0.01, type=float, help='Binwidth for DOS histograms. Unit: eV.')
parser.add_argument('--trace', default=0, type=int, help='Numbers of HOMO states to trace.')
parser.add_argument('--reduce', action='store_true', help='Joining all orbitals into one occupation.')

class Cp2kPdosFrame(object):
	_iterstep = None
	_fermi = None
	_data = None

	def __init__(self, lines):
		self._labels = []
		try:
			self._iterstep = int(lines[0].split('i = ')[1].split(',')[0])
		except:
			self._iterstep = int(lines[0].split()[-3])
		try:
			self._fermi = float(lines[0].split()[-2])
		except:
			self._fermi = None

		self._labels = lines[1].split()[5:]
		self._data = np.zeros((len(self._labels)+3, len(lines)-2))
		for idx, line in enumerate(lines[2:]):
			self._data[:, idx] = map(float, line.split())

	def join(self, other):
		for idx, orbital in enumerate(other._labels):
			try:
				this_idx = self._labels.index(orbital)
			except:
				this_idx = len(self._labels)
				self._labels.append(orbital)
				data = np.zeros((self._data.shape[0]+1, self._data.shape[1]))
				data[:-1, :] = self._data
				self._data = data
			self._data[this_idx+3, :] += other._data[idx+3, :]

	def _get_indices(self, orbital):
		if orbital == '*':
			return np.array(range(len(self._labels)))

		if orbital in self._labels:
			return np.array([self._labels.index(orbital)])

		if len(orbital) == 1:
			idxs = []
			for idx, label in enumerate(self._labels):
				if label.startswith(orbital):
					idxs.append(idx)
			if len(idxs) == 0:
				raise ValueError()
			return np.array(idxs)

		raise ValueError()

	def _count_orbitals(self):
		return self._data.shape[1]

	def invert(self):
		self._data[1:, :] *= -1		

	def center(self, orbital, cutoff, relative_fermi):
		try:
			orbital, selection = orbital.split('.')
		except:
			selection = 'all'
		try:
			idxs = self._get_indices(orbital)
		except ValueError:
			raise ValueError('Orbital %s not in iteration step %d' % (orbital, self._iterstep))

		# filter entries that are to be ignored based on occupation and cutoff value
		validfilter = np.ones(self._count_orbitals(), dtype=bool)
		if 'all'.startswith(selection):
			# include both occupied and unoccupied orbitals
			pass
		if 'valence'.startswith(selection):
			# include valence band only
			validfilter = self._data[2] == 1
		if 'conduction'.startswith(selection):
			# include conduction band only
			validfilter = self._data[2] == 0
		validfilter *= self._data[1] > cutoff

		weights = np.sum(self._data[idxs+3, :]*validfilter, axis=1)
		vals = np.sum(self._data[idxs+3, :]*validfilter*self._data[1], axis=1)
		return (np.sum(vals) / np.sum(weights)) - relative_fermi*self._fermi
	
	def histogram(self, orbital, bin, cutoff, relative_fermi):
		try:
			orbital, selection = orbital.split('.')
		except:
			selection = 'all'
		try:
			idxs = self._get_indices(orbital)
		except ValueError:
			raise ValueError('Orbital %s not in iteration step %d' % (orbital, self._iterstep))

		
		validfilter = self._data[1] > cutoff
		energies = self._data[1] - relative_fermi * self._fermi
		weights = np.sum(self._data[idxs+3, :], axis=0)
		mval = max(cutoff, energies.min())
		weights *= validfilter

		if isinstance(bin, float):
			bins = np.linspace(mval, energies.max(), (energies.max()-mval)/bin)
		else:
			bins = np.array(bin)

		hs, bins = np.histogram(energies, bins, (bins.min(), bins.max()), weights=weights, density=True)
		return hs, bins

	def iterstep(self):
		return self._iterstep

	def get_trace(self, nhomo, relative_fermi, relative_homo):
		rows = np.where(self._data[2, :] == 1)[0][-nhomo:]
		sel = self._data[:, rows]
		energies = sel[1] - relative_fermi*self._fermi - relative_homo*sel[1][-1]
		occupations = np.sum(sel[3:], axis=0)
		return energies[::-1], occupations[::-1]

	def reduce(self, relative_fermi):
		return self._data[1]-relative_fermi*self._fermi, np.sum(self._data[3:], axis=0)

class Cp2kPdosFile(object):
	def __init__(self, fh, limit=0, skip=0):
		self._frames = []
		t_lines = []

		if fh.name[-3:] == '.gz':
			fh = gzip.GzipFile(fileobj=fh)

		read_frames = 0
		for line in fh:
			if line.startswith('# ') and len(t_lines) > 1:
				read_frames += 1
				if read_frames > skip:
					self._frames.append(Cp2kPdosFrame(t_lines))
				t_lines = []

			if read_frames == limit and limit != 0:
				break

			t_lines.append(line)
		if (read_frames < limit or limit == -1) and read_frames > skip:
			self._frames.append(Cp2kPdosFrame(t_lines))
		fh.close()

	def __len__(self):
		return len(self._frames)

	def list_orbitals(self):
		if len(self) == 0:
			raise ValueError('Empty file.')

		return self._frames[0]._labels

	def center(self, orbitals, cutoff, relative_fermi):
		orbitals = orbitals.split()
		pos = np.zeros((len(self._frames), len(orbitals)))
		for odx, orbital in enumerate(orbitals):
			for idx, frame in enumerate(self._frames):
				pos[idx, odx] = frame.center(orbital, cutoff, relative_fermi)
		return [_.iterstep() for _ in self._frames], pos * 27.21138505

	def join(self, other):
		for idx in range(len(other._frames)):
			self._frames[idx].join(other._frames[idx])
	
	def dos(self, bin, cutoff, relative_fermi, smear):
		orbitals = '*'
		hs, bins = self._frames[0].histogram('*', bin, cutoff, relative_fermi)
		for frame in self._frames[1:]:
			t, _ = frame.histogram('*', bins, cutoff, relative_fermi)
			hs += t
		hs /= len(self)
		xs = (bins[1:]+bins[:-1])/2
		xs2 = xs - np.average(xs)
		gs = np.exp(-xs2**2/(2*smear**2))/(smear*np.sqrt(2*np.pi))
		smeared = np.convolve(hs, gs, mode='same')
		return xs, hs, smeared

	def projectout(self, other):
		for idx in range(len(other._frames)):
			other._frames[idx].invert()
		self.join(other)

	def take_fermi_from(self, other):
		if len(self._frames) != len(other._frames):
			raise TypeError('Incompatible files')

		for idx in range(len(self._frames)):
			self._frames[idx]._fermi = other._frames[idx]._fermi

	def trace(self, nhomo, relative_fermi, relative_homo):
		if relative_homo and relative_fermi:
			raise ValueError('Only one reference value allowed: got Fermi and HOMO')
		results = np.zeros((len(self), 1+nhomo*2))

		for frame in range(len(self)):
			results[frame, 0] = frame
			energies, occupations = self._frames[frame].get_trace(nhomo, relative_fermi, relative_homo)
			results[frame, 1::2] = energies
			results[frame, 2::2] = occupations
		return results

	def reduce(self, relative_fermi):
		if len(self) > 1:
			print '# WARNING: ONLY TAKING FIRST FRAME INTO ACCCOUNT'
		return self._frames[0].reduce(relative_fermi)

def main():
	args = parser.parse_args()
	names = [args.file.name]
	pdos = Cp2kPdosFile(args.file, limit=args.limit)
	# fermi level from seperate file?
	if args.relative_fermi and args.fermi_file != None:
		fermi = Cp2kPdosFile(args.fermi_file, limit=args.limit, skip=args.skip)
		pdos.take_fermi_from(fermi)

	# merge / split in the beginning
	if args.add is not None:
		for pdosfile in args.add:
			tpdos = Cp2kPdosFile(pdosfile, limit=args.limit, skip=args.skip)
			pdos.join(tpdos)
			names.append(pdosfile.name)
	if args.subtract is not None:
		for pdosfile in args.subtract:
			tpdos = Cp2kPdosFile(pdosfile, limit=args.limit, skip=args.skip)
			pdos.projectout(tpdos)
			names.append(pdosfile.name)

	# print metadata
	print '# Read file(s) %s with %d iteration steps' % (', '.join(names), len(pdos))
	print '# Found orbitals: %s' % ' '.join(pdos.list_orbitals())
	print '# Selected centers: %s' % args.orbital_center
	print '# Orbitals ignored below: %f Hartree' % args.orbital_cutoff
	print '# Output unit: eV'
	if args.orbital_center is not None:
		itersteps, data = pdos.center(args.orbital_center, args.orbital_cutoff, args.relative_fermi)
		for i, e in zip(itersteps, data):
			print i, 
			for t in e:
				print t,
			print
		if args.visualise:
			import matplotlib.pyplot as plt 
			scount, ocount = data.shape
			for idx in range(ocount):
				plt.plot(range(scount), data[:, idx])
			plt.show()
	if args.dos:
		xs, hs, gs = pdos.dos(args.binwidth/27.21138505, args.dos_cutoff, args.relative_fermi, args.dos_smear/27.21138505)		
		xs *= 27.21138505
		for x, h, g in zip(xs, hs, gs):
			print x, h, g
		if args.visualise:
			import matplotlib.pyplot as plt
			plt.plot(xs, -hs/hs.max())
			plt.plot(xs, gs/gs.max())
			plt.show()
	if args.trace != 0:
		results = pdos.trace(args.trace, args.relative_fermi, args.relative_homo)
		results[:, 1::2] *= 27.21138505
		print '# Frame Energy(HOMO 1, highest HOMO), Occupation(HOMO 1), Energy(HOMO 2), ...'
		for frame in results:
			for e in frame:
				print e,
			print 
	if args.reduce != 0:
		energies, res = pdos.reduce(args.relative_fermi)
		energies *= 27.21138505
		for i, e in enumerate(zip(energies, res)):
			e, v = e
			print i, e, v

if __name__ == '__main__':
	main()

