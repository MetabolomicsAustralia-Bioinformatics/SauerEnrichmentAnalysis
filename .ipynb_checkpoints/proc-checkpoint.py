#from main_script import read_mhunter_csv, calculate_labelling, read_atomic_composition
from SauerFunction import *

from SauerClass import Record, Labelling


# Read all the relevant files
records = read_mhunter_csv('sample_files/sample_input_data.csv')
atomic_composition, N_dict = read_atomic_composition('sample_files/sample_input_formulae.csv')
fp = open('sample_files/out_v2.csv', 'w')


labelling_list = []
max_results_length = 0

for record in records:
    results_dict = calculate_labelling(record, N_dict, atomic_composition)
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
print("Written out results.")