#!/usr/env python3
import SauerFunction as sf
from SauerClass import Record, Labelling
import argparse
import sys

VERSION = '0.1'
AUTHORS = 'Sean O\'Callaghan, Michael Leeming, Don Teng'


parser = argparse.ArgumentParser(description='Do enrichment analysis.')
parser.add_argument('-1','--input_data', help='Input data file name')
parser.add_argument('-2','--input_formulae', help='Input formulae file name')
parser.add_argument('-o','--output_filename', help='Name of output file')
parser.add_argument('--version', help='Print version and exit', action='store_true')

args = parser.parse_args()

if args.version:
    print('proc.py %s' % VERSION, file=sys.stderr)
    exit(0)

records = sf.read_mhunter_csv(args.input_data)
atomic_composition, N_dict = sf.read_atomic_composition(args.input_formulae)
fp = open(args.output_filename, 'w')


labelling_list = []
max_results_length = 0

for record in records:
    results_dict = sf.calculate_labelling(record, N_dict, atomic_composition)
    for key, value in results_dict.items():
        if len(value) > max_results_length:
            max_results_length = len(value)
    labelling_list.append(Labelling(record.get_name(), results_dict))

fp.write("Species, Labelling Source, Sample Name, Labelling %,")
for i in range(max_results_length-1):
    fp.write("m" + str(i) + ",")
fp.write("\n")

for label in labelling_list:
    species = label.get_species()
    label_dict = label.get_label_dict()
    names = label_dict.keys()
    names = sorted(names)
    i = 0
    for name in names:
        for key, value in label_dict.items():
            if name == key:
                fp.write(species + ',')
                if len(key.split(',')) == 2:
                    fp.write(key.split(',')[1].strip('"') + ',')
                    fp.write(key.split(',')[0].strip('"') + ',')
                else:
                    fp.write(' ,')
                    fp.write(key + ',')
                for val in value:
                    fp.write(str(val) + ', ')
                fp.write('\n')
print("Process complete; wrote out results.")

