# Convert aracne2 .adj files to tab separated values suitable .sif files for cytoscape
# Usage: $ python adj2cytoscape.py file.adj > file.sif

import argparse

parser = argparse.ArgumentParser(description='Convert aracne2 adjacency files to tab separated values suitable for cytoscape.')
parser.add_argument('adj', type=argparse.FileType('r'), nargs="+", help="One or more ADJ files." )
args    = parser.parse_args()

for f in args.adj:
    for l in f.readlines():
        if not l.startswith('>'):
            cols   = iter(l.strip().split('\t'))
            source = cols.next()
            for target, mi in ((item, cols.next()) for item in cols):
                print "%s %s %s" % (source, mi, target)

