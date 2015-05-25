#!/usr/bin/env python

import argparse
import euston.io as io
import euston.helper as helper

parser = argparse.ArgumentParser(description='Performance analysis of CP2K output files.')
parser.add_argument('input', type=str, help='CP2K log file')

def main(parser):
    args = parser.parse_args()

    cp2k = io.Cp2kLog(args.input)
    cores = cp2k.get_num_cores()
    spc = cp2k.get_num_spc()
    time = cp2k.get_time_spc()
    mds = cp2k.get_num_md()

    print 'Number of cores:     {0:>25}'.format(cores)
    print 'Number of SCF steps: {0:>25}'.format(spc)
    print 'SCF time:            {0:>25} ({1}s)'.format(helper.human_readable_time(time), time)
    print 'MD steps:            {0:>25}'.format(mds)
    print
    print 'SCF steps / MD step: {0:>25}'.format(spc*1.0/mds)
    timepermd = time / mds
    print 'SCF time / MD step:  {0:>25} ({1}s)'.format(helper.human_readable_time(timepermd), timepermd)
    timeperscf = time / spc
    print 'SCF time / SCF step: {0:>25} ({1}s)'.format(helper.human_readable_time(timeperscf), timeperscf)


if __name__ == '__main__':
    main(parser)
