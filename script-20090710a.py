"""
Remove some columns from the data file.

The first column is not removed.
The number of non-first columns to keep is specified (e.g. 40)
The kept columns are chosen randomly
"""

import random
import sys

class UsageError(Exception): pass

def main():
    usage = 'Usage: python ' + sys.argv[0] + ' <ncols> <infile> <outfile>'
    try:
        # read the command line arguments
        if len(sys.argv) != 4:
            print 'expected three arguments'
            raise UsageError(usage)
        scriptname, ncols_string, inpath, outpath = sys.argv
        # get the number of columns to keep
        try:
            ncols = int(ncols_string)
        except:
            raise UsageError(usage)
        # read the input file
        fin = open(inpath)
        lines = fin.readlines()
        fin.close()
        # skip empty lines
        lines = [line for line in lines if line]
        # get rows of elements
        rows = [[x.strip() for x in line.split(',')] for line in lines]
        # get the columns of elements
        columns = zip(*rows)
        nobsdatacols = len(columns) - 1
        if nobsdatacols < ncols:
            print 'the input file has only', nobsdatacols, 'columns of data'
            raise UsageError(usage)
        # choose a subset of the columns
        columns = [columns[0]] + random.sample(columns[1:], ncols)
        # convert back to rows
        rows = zip(*columns)
        # write the output file
        fout = open(outpath, 'w')
        for row in rows:
            print >> fout, ','.join(row)
        fout.close()
    except UsageError, e:
        print e

if __name__ == '__main__':
    main()
