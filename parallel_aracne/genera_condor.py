import argparse
import os, stat
from jinja2 import Environment, FileSystemLoader

parser = argparse.ArgumentParser(description='Generates condor submit file for aracne runs.')
parser.add_argument('--path_to_aracne2',  type=argparse.FileType('r'), required=True, help='path to aracne2 binary')
parser.add_argument('--expfile',  type=argparse.FileType('r'), required=True, help='expression file')
parser.add_argument('--probes',  type=argparse.FileType('r'), required=True, help='probes, one in every line' )
parser.add_argument('--run_id',  required=True, help="name of condor run" )
parser.add_argument('--outdir',  required=True, help="outdir for adj matrices" )
parser.add_argument('--p',  required=True, help="P-value: e.g. 1e-7" )

args    = parser.parse_args()

expfile = args.expfile.name
p       = args.p
outdir  = args.outdir


# make sane affy ids
probes = []
for id in args.probes.readlines():
    probes.append( id.strip() )


# create exec dir
if not os.path.exists(outdir):
    os.makedirs(outdir)


# create aracne.sh
aracne_path = os.path.dirname(os.path.realpath(args.path_to_aracne2.name))
aracne_sh = """#!/bin/bash
cd {aracne_path}
./aracne2 $@"""
with open(os.path.join(outdir,'aracne.sh'), 'w') as f:
    f.write(aracne_sh.format(aracne_path=aracne_path))
    os.chmod(f.name, stat.S_IREAD | stat.S_IWRITE | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH )





#
# create condor script
#
scriptname = "%s/%s.condor" % (outdir, args.run_id)

# use same dir as this file's as environment
env = Environment(loader=FileSystemLoader(os.path.dirname(os.path.realpath(__file__))))
template = env.get_template('condor_aracne.tt')

with open(scriptname, 'w') as f:
    f.write( template.render( expfile = expfile,
                              probes = probes,
                              p       = p,
                              outdir  = outdir,
                              run_id  = args.run_id ) )
