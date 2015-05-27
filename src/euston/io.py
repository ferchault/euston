#!/usr/bin/env python

# system modules
import abc
import gzip
import math
import itertools
import re

# third-party modules
import numpy as np
try:
    import MDAnalysis as mda
    HAS_MDA = True
except:
    HAS_MDA = False

# custom modules
import geometry as geo

BOHR2ANGSTROM = 1/0.529177210

def require_parsed(f):
    def _require_parsed(self, *args, **kwargs):
        if not self._parsed:
            raise AssertionError('Function only available after completion of parsing the input file.')
        return f(self, *args, **kwargs)
    return _require_parsed

def require_loaded(f):
    def _require_loaded(self, *args, **kwargs):
        if not self._loaded:
            raise AssertionError('Function only available after loading an input file or similar input data.')
        return f(self, *args, **kwargs)
    return _require_loaded

def write_lines(filename, lines, terminate='\n'):
    try:
        fh = open(filename, 'w')
    except:
        raise ValueError('Unable to open file for writing.')

    for line in lines:
        fh.write('%s%s' % (line, terminate))

    fh.close()

class HoldsCoordinates(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_coordinates(self):
        pass

class HoldsUnitcell(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_h_matrix(self):
        pass

def anyopen(filename):
    rawname = filename
    gzip = '.gz .gzip'.split()
    for ext in gzip:
        if rawname.endswith(ext):
            rawname = rawname[:-len(ext)]

    if rawname.endswith('.xyz'):
        return XYZ(filename=filename)

    if rawname.endswith('.dcd'):
        return DCD(filename=filename)

class FileIO(object):
    """Abstract base class for all file objects with support for gzipped input files.."""
    
    #: Whether the instance has been populated with data already
    _loaded = False
    #: File handle for input file, if applicable
    _fh = None
    #: Whether the input has been parsed already
    _parsed = False

    def __init__(self, filename=None, filehandle=None):
        """Prepares reading input files.

        :param filename: Optional input filename.
        """
        if filename is not None and filehandle is not None:
            raise ValueError('Only one argument (either filename or filehandle) is allowed.')
        if filename is None and filehandle is None:
            #: Whether the instance has been populated with data already
            self._loaded = False
            #: File handle for input file, if applicable
            self._fh = None
        elif filehandle is not None:
            self._fh = filehandle
            self._loaded = True
        else:
            self._loaded = True
            if filename[-3:] == '.gz' or filename[-5:] == '.gzip':
                self._fh = gzip.open(filename, 'rb')
            else:
                self._fh = open(filename, 'r')

        #: Whether the input has been parsed already
        self._parsed = False
        self._parse()

    def _parse(self):
        """Finalises file content parsing."""
        self._parsed = True

class DCD(HoldsCoordinates, HoldsUnitcell, FileIO):
    @require_loaded
    def _parse(self):
        if not HAS_MDA:
            raise TypeError('MDAnalysis required in order to load a DCD file.')
        reader = mda.coordinates.DCD.DCDReader(self._fh.name)
        self._hmat = geo.abc_to_hmatrix(*reader.ts.dimensions, degrees=True)
        self._coordinates = np.zeros((reader.ts.numatoms, 3))
        self._coordinates[:, 0] = reader.ts._x
        self._coordinates[:, 1] = reader.ts._y
        self._coordinates[:, 2] = reader.ts._z

        super(DCD, self)._parse()

    @require_loaded
    @require_parsed
    def get_coordinates(self):
        return self._coordinates

    @require_loaded
    @require_parsed
    def get_h_matrix(self):
        return self._hmat

class XYZ(HoldsCoordinates, FileIO):
    @require_loaded
    def count_atoms(self):
        return self._coordinates.shape[0]

    def set_data(self, labels, coord):
        if len(labels) != coord.shape[0]:
            raise ValueError('Mismatching lengths for labels and coordinates.')

        self._coordinates = coord
        self._labels = labels

    def to_string(self):
        lines = []
        lines.append('%d' % len(self._labels))
        lines.append(self._comment)
        for label, atom in zip(self._labels, self._coordinates):
            lines.append('%s %f %f %f' % (label, atom[0], atom[1], atom[2]))
        return lines

    def _parse(self):
        if self._fh is None:
            self._comment = ''
            return
        lines = self._fh.readlines()
        try:
            num_atoms = int(lines[0].strip())
        except:
            raise ValueError('Invalid atom count specified.')
        if len(lines) < num_atoms+2:
            raise ValueError('XYZ file contains less atoms than specified.')
        if len(lines) > num_atoms+2:
            raise NotImplementedError('Only single frame XYZ files are supported.')

        self._coordinates = np.zeros((num_atoms, 3))
        self._labels = []
        self._comment = lines[1].strip()
        for line in lines[2:]:
            parts = line.strip().split()
            if len(parts) != 4:
                raise ValueError('Invalid atom line: %s' % line.strip())
            self._labels.append(parts[0])
            try:
                coor = map(float, parts[1:])
            except:
                raise ValueError('Invalid coordinates: %s' % line.strip())
            idx = len(self._labels)-1
            self._coordinates[idx, :] = np.array(coor)

        super(XYZ, self)._parse()

    def get_coordinates(self):
        return self._coordinates

    def get_labels(self):
        return self._labels

class Cp2kLog(FileIO):
    def get_values_matching(self, regex, count=None, line_numbers=False, transform=(lambda x: x, )):
        result = []
        for no, line in enumerate(self._lines):
            found = re.findall(regex, line)
            if len(found) != 0:
                tres = [[func(val) for func, val in zip(transform, found)]]
                if line_numbers:
                    tres.append(no)
                result.append(tres)
                if count is not None and len(result) == count:
                    break
        return result

    @require_parsed
    def get_num_cores(self):
        if self._num_cores is not None:
            return self._num_cores

        pattern = " GLOBAL\| Total number of message passing processes[ ]*(\d*)"
        vals = self.get_values_matching(pattern, count=1, transform=(int, ))
        self._num_cores = vals[0][0][0]
        return self._num_cores

    @require_parsed
    def get_num_spc(self):
        if self._num_spc is not None:
            return self._num_spc

        pattern = "SCF run converged in[ ]*(\d*) steps"
        vals = self.get_values_matching(pattern, transform=(int, ))
        self._num_spc = 0
        for val in vals:
            self._num_spc += val[0][0]
        return self._num_spc

    @require_parsed
    def get_time_spc(self):
        if self._time_spc is not None:
            return self._time_spc

        pattern = "CPU TIME \[s\]                 =[ ]*(\d*\.\d*).*"
        vals = self.get_values_matching(pattern, transform=(float, ))
        self._time_spc = 0
        for val in vals:
            self._time_spc += val[0][0]
        return self._time_spc

    @require_parsed
    def get_num_md(self):
        if self._num_md is not None:
            return self._num_md

        pattern = "STEP NUMBER"
        vals = self.get_values_matching(pattern)
        self._num_md = len(vals)
        return self._num_md

    def _parse(self):
        self._lines = self._fh.readlines()
        runs = self.get_values_matching('.* PROGRAM STARTED AT .*', line_numbers=True)
        if len(runs) < 0:
            raise ValueError('No CP2K run detected. Damaged log file?')
        if len(runs) > 1:
            print 'Warning: Only keeping last run.'
            self._lines = self._lines[runs[-1][1]:]

        self._num_cores = None
        self._time_spc = None
        self._num_spc = None
        self._num_md = None
        super(Cp2kLog, self)._parse()

class Cp2kInput(FileIO):
    def set_lines(self, lines):
        self._lines = lines
        self._prepare()

    @require_loaded
    def to_string(self, close_sections=True, indent='  ', keep_comments=True, indent_comments=True, keep_empty=False):
        lines = []
        sections = []
        line_numbers = []
        for no, line in enumerate(self._lines):
            line = line.strip()
            if len(line) == 0:
                if keep_empty:
                    lines.append('')
                continue
            if (line.startswith('#') or line.startswith('!')):
                if not keep_comments:
                    continue
                if indent_comments:
                    lines.append('%s%s' % (indent*len(sections), line))
                else:
                    lines.append(line)
                continue
            if line.startswith('&END'):
                ending = line[4:].strip()
                if ending != '' and (len(sections) == 0 or ending != sections[-1]):
                    raise ValueError('Trying to end the %s section from line %d with %s in line %d' % (sections[-1], line_numbers[-1]+1, ending, no+1))
                if close_sections:
                    lines.append('%s&END %s' % (indent*(len(sections)-1), sections[-1]))
                else:
                    lines.append('%s&END' % (indent*(len(sections)-1)))
                lines.append('')
                sections = sections[:-1]
                line_numbers = line_numbers[:-1]
                continue
            if line.startswith('&'):
                if len(lines) != 0:
                    lines.append('')
                lines.append('%s%s' % (indent*len(sections), line))
                sections.append((line.split()[0][1:]).strip())
                line_numbers.append(no)
                continue
            lines.append('%s%s' % (indent*len(sections), line))

        retlines = [lines[0]]
        for line in lines[1:]:
            if line == '':
                if retlines[-1] == '':
                    continue
            if line.strip().startswith('&END ') and retlines[-1] == '' and len(retlines) > 1 and retlines[-2].strip().startswith('&END '):
                retlines[-1] = line
                continue
            if line.strip().startswith('&') and not line.strip().startswith('&END ') and retlines[-1] == '' and len(retlines) > 1 and retlines[-2].strip().startswith('&') and not retlines[-2].strip().startswith('&END '):
                retlines[-1] = line
                continue
            retlines.append(line)

        if retlines[-1] == '':
            retlines.pop()
        return retlines

    @require_loaded
    def get_keyword_checked(self, keyword, conversion=lambda x: x):
        try:
            res = self.get_path(keyword)
        except ValueError as e:
            #print str(e)
            return None
        if res is None:
            #print 'No such key: %s' % keyword
            return None
        try:
            val = map(conversion, res.split())
        except:
            #print 'Invalid %s entry.' % keyword
            return None
        return val

    def boolean(self, val, default):
        if val is None:
            return default
        if val in ('T', '.T.', 'TRUE'):
            return True
        return False

    @require_loaded
    def get_path(self, path):
        elements = path.split(' / ')
        sections = []
        line_numbers = []
        collect = []
        for no, line in enumerate(self._lines):
            if len(line) == 0 or line.startswith('#') or line.startswith('!'):
                continue
            if line.startswith('&END'):
                ending = line[4:].strip()
                if ending != '' and ending != sections[-1]:
                    raise ValueError('Trying to end the %s section from line %d with %s in line %d' % (sections[-1], line_numbers[-1]+1, ending, no+1))
                sections = sections[:-1]
                line_numbers = line_numbers[:-1]
                continue
            if line.startswith('&'):
                sections.append(line.split()[0][1:])
                line_numbers.append(no)
                continue

            if sections == elements[:-1]:
                if elements[-1] == '*':
                    collect.append(line)
                    continue
                parts = line.split()
                if parts[0] == elements[-1]:
                    return ' '.join(parts[1:])
        if elements[-1] == '*':
            return collect
        return None

    def _prepare(self):
        self._lines = [_.strip() for _ in self._lines]
        self._loaded = True

    def _parse(self):
        if self._fh is not None:
            self._lines = self._fh.readlines()
            self._prepare()

        # finalise parsing
        super(Cp2kInput, self)._parse()

    @require_loaded
    def get_cell_vectors(self):
        a, b, c = (None, None, None)
        alpha, beta, gamma = (None, None, None)

        # check whether ABC / ALPHA_BETA_GAMMA is set
        retval1 = self.get_keyword_checked('FORCE_EVAL / SUBSYS / CELL / ABC', float)
        retval2 = self.get_keyword_checked('FORCE_EVAL / SUBSYS / CELL / ALPHA_BETA_GAMMA', float)
        if retval1 is not None:
            if retval2 is None:
                retval2 = [90, 90, 90]
            a, b, c = retval1
            alpha, beta, gamma = np.radians(np.array(retval2))

            # convert to vectors
            hmat = geo.abc_to_hmatrix(a, b, c, alpha, beta, gamma, degrees=False)
            a = hmat[:, 0]
            b = hmat[:, 1]
            c = hmat[:, 2]

        # check whether A, B, C is set
        if a is None:
            retval1 = self.get_keyword_checked('FORCE_EVAL / SUBSYS / CELL / A', float)
            retval2 = self.get_keyword_checked('FORCE_EVAL / SUBSYS / CELL / B', float)
            retval3 = self.get_keyword_checked('FORCE_EVAL / SUBSYS / CELL / C', float)
            if retval1 is not None and retval2 is not None and retval3 is not None:
                a = np.array(retval1)
                b = np.array(retval2)
                c = np.array(retval3)

        if a is None:
            print 'No supported cell information found.'
            return None

        return (a, b, c)


class CubeFile(HoldsUnitcell, FileIO):
    def count_atoms(self):
        return self._natoms

    def count_voxels(self):
        return abs(reduce(lambda x, y: x*y, self._nvoxel))

    @require_loaded
    @require_parsed
    def get_h_matrix(self):
        h = np.copy(self._vectors).transpose()
        for i in range(3):
            h[:, i] *= self._nvoxel[i]

        return h

    @require_loaded
    @require_parsed
    def get_val(self, x, y, z):
        return self._data[x, y, z]

    @require_loaded
    @require_parsed
    def get_coordinates(self):
        return np.copy(self._coordinates)

    @require_loaded
    @require_parsed
    def set_coordinates(self, coord):
        try:
            this_shape = coord.shape
        except:
            # not a numpy array
            coord = np.array(coord).reshape((-1, 3))
            this_shape = coord.shape
        if this_shape != self._coordinates.shape:
            raise ValueError('Changing coordinate shape not implemented.')

        self._coordinates = np.copy(coord)

    @require_loaded
    @require_parsed
    def get_xlen(self):
        return self._nvoxel[0]

    @require_loaded
    @require_parsed
    def get_ylen(self):
        return self._nvoxel[1]

    @require_loaded
    @require_parsed
    def get_zlen(self):
        return self._nvoxel[2]

    @require_loaded
    @require_parsed
    def get_voxel_pos(self, x, y, z, centered=False):
        if centered:
            x += .5
            y += .5
            z += .5
        return x*self._vectors[0, :] + y*self._vectors[1, :] + z * self._vectors[2, :]

    @require_loaded
    @require_parsed
    def get_voxel_volume(self):
        return geo.cell_volume(np.copy(self._vectors).transpose())

    @require_loaded
    @require_parsed
    def get_projection(self, axis_index, absolute):
        other_axes = set(range(3)) - set([axis_index])
        if absolute:
            t = np.sum(np.absolute(self._data), axis=max(other_axes))
        else:
            t = np.sum(self._data, axis=max(other_axes))
        return np.sum(t, axis=min(other_axes))

    @require_loaded
    def _parse(self):
        # default values
        self._origin = np.array([0, 0, 0])
        self._natoms = 0
        self._vectors = np.zeros((3, 3))
        self._nvoxel = np.array([0, 0, 0])
        self._atomic_numbers = []
        self._coordinates = None
        self._data = None
        self._header = []

        # file length
        lines = []
        while len(lines) < 2+1+3:
            try:
                lines.append(next(self._fh).strip())
            except:
                raise ValueError('File too short.')

        # headers
        self._header = lines[:2]

        # atom count and origin
        parts = lines[2].split()
        try:
            parts[0] = int(parts[0])
            parts[1:] = map(float, parts[1:])
        except:
            parts = [None]
        if len(parts) != 4:
            raise ValueError('Invalid definition of atom count and origin.')
        self._origin = np.array(parts[1:])
        self._natoms = parts[0]

        # number of voxels and voxel vectors
        parts = map(str.split, lines[3:6])
        try:
            for idx, part in enumerate(parts):
                if len(part) != 4:
                    raise ValueError()
                self._nvoxel[idx] = int(part[0])
                self._vectors[idx, :] = map(float, part[1:])

                if self._nvoxel[idx] > 0:
                    self._vectors[idx, :] /= BOHR2ANGSTROM
        except:
            raise ValueError('Invalid voxel size definition.')

        # positions
        adata = np.genfromtxt(itertools.islice(self._fh, self._natoms))
        if adata.ndim == 1:
            self._atomic_numbers = [adata[0]]
            self._coordinates = np.array([adata[2:]])
        else:
            self._atomic_numbers = adata[:, 0]
            self._coordinates = adata[:, 2:]
        for idx, val in enumerate(self._nvoxel):
            if val > 0:
                self._coordinates[:, idx] /= BOHR2ANGSTROM

        # voxel
        self._data = np.zeros(self.count_voxels())
        count = 0
        try:
            for line in self._fh: #lines[(6 + self._natoms):]:
                for number in map(float, line.split()):
                    self._data[count] = number
                    count += 1
        except:
            raise
        if count != len(self._data):
            raise ValueError('Truncated voxel data.')

        # rescaling values in case the axes have been rescaled
        self._data /= BOHR2ANGSTROM**(np.sum(self._nvoxel > 0))
        self._originalnvoxel = np.copy(self._nvoxel)
        self._nvoxel = np.abs(self._nvoxel)
        self._data = np.reshape(self._data, tuple(self._nvoxel))

        # finalise parsing
        super(CubeFile, self)._parse()

    @require_loaded
    @require_parsed
    def to_string(self):
        lines = self._header
        lines.append('%d %e %e %e' % (self._natoms, self._origin[0], self._origin[1], self._origin[2]))
        scale = np.ones(3)
        for i in range(3):
            if self._originalnvoxel[i] > 0:
                scale[i] = BOHR2ANGSTROM
            lines.append('%d %f %f %f' % (self._originalnvoxel[i], self._vectors[i, 0]*scale[i], self._vectors[i, 1]*scale[i], self._vectors[i, 2]*scale[i]))
        for i in range(len(self._atomic_numbers)):
            coord = scale*self._coordinates[i]
            lines.append('%d 0 %f %f %f' % (self._atomic_numbers[i], coord[0], coord[1], coord[2]))
        scale = np.prod(scale)
        segment = []
        for x in range(self._nvoxel[0]):
            for y in range(self._nvoxel[1]):
                for z in range(self._nvoxel[2]):
                   segment.append(self._data[x, y, z]*scale)
                   if len(segment) == 5:
                       lines.append('%e %e %e %e %e' % (segment[0], segment[1], segment[2], segment[3], segment[4]))
                       segment = []
        lines.append(' '.join(map(lambda _: '%e' % _, segment)))

        return lines
