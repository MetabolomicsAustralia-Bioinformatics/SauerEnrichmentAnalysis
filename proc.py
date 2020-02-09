#from main_script import read_mhunter_csv, calculate_labelling, read_atomic_composition
import SauerFunction as sf
from SauerClass import Record, Labelling
import argparse



# Read all the relevant files
#records = sf.read_mhunter_csv('sample_files/sample_input_data.csv')
#atomic_composition, N_dict = sf.read_atomic_composition('sample_files/sample_input_formulae.csv')
#fp = open('sample_files/out_temp.csv', 'w')


parser = argparse.ArgumentParser(description='Do enrichment analysis.')
parser.add_argument('-1','--input_data', help='Input data file name', required=True)
parser.add_argument('-2','--input_formulae', help='Input formulae file name', required=True)
parser.add_argument('-o','--output_filename', help='Name of output file', required=True)
args = vars(parser.parse_args())

records = sf.read_mhunter_csv(args["input_data"])
atomic_composition, N_dict = sf.read_atomic_composition(args["input_formulae"])
fp = open(args["output_filename"], 'w')


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

