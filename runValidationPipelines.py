#!/cm/shared/apps/python/3.6.5/bin/python
import argparse
import shutil
import os
import subprocess

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

class CustomFormatter(argparse.RawDescriptionHelpFormatter):
    pass

epilog ="""
    Steps performed by this script:
    1. check to see if BIP outputs exist, or else schedule bip runs
    2. Run scripts over the different samples and different bip versions
    3. create an summary tables for the variants of interest in the study
    4. create a reportable table comparing the values between the two versions of bip
    """

parser = argparse.ArgumentParser(description="BIP accuracy comparison runner",epilog=epilog, formatter_class=CustomFormatter)
parser.add_argument("-o", "--output_dir",
                    action="store",
                    dest="output_dir",
                    required=True,
                    help=("directory to store the tables"))
parser.add_argument("-d", "--_data_tables_dir",
                    action="store",
                    dest="_data_tables_dir",
                    required=True,
                    help=("directory of the sample and variant tables for the study"))
parser.add_argument("-v1", "--version1",
                    action="store",
                    dest="version1",
                    required=True,
                    help=("specify the version of the BIP used"))
parser.add_argument("-v2", "--version2",
                    action="store",
                    dest="version2",
                    required=True,
                    help=("specify the version of the BIP used"))

args = parser.parse_args()

if not os.path.isdir(args._data_tables_dir):
    raise RuntimeError("Nonexistant or invalid input directory '{}'".format(args.in_dir))

abs_dir = subprocess.Popen(["pwd"], stdout=subprocess.PIPE, shell=True).communicate()[0].decode('utf-8').replace('\n','/')

out_dir = args.output_dir
data_dir = args._data_tables_dir

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

bip_out_dir = "/ghds/ivd/flowcentral/"

version1, version2 = args.version1, args.version2

def get_version(version):
    if version in ["3.5", "3.5.0", "3_5", "3_5_0"]:
        return "3.5.0"
    if version in ["3.5.1", "3_5_1"]:
        return "3.5.1"
    if version in ["3.5.2", "3_5_2"]:
        return "3.5.2"
    if version in ["3.5.3", "3_5_3"]:
        return "3.5.3"
    else:
        return version
    
version1, version2 = get_version(version1), get_version(version2)

print('####################')
print('##  compare accuracy')
print('####################')

accuracy_data_dir = data_dir + "accuracy/"
accuracy_out_dir = out_dir + "accuracy_comparsion/"

if not os.path.isdir(accuracy_out_dir):
    os.mkdir(accuracy_out_dir)

os.system("python3 BIP_validation_pipelines_accuracy/compare_accuracy_values.py -o " + accuracy_out_dir + " -d1 " + accuracy_data_dir + version1+"/" + " -v1 " + version1 + " -d2 " + accuracy_data_dir + version2+"/" + " -v2 " + version2)
#print("python3 BIP_validation_pipelines_accuracy/compare_accuracy_values.py -o " + accuracy_out_dir + " -d1 " + accuracy_data_dir + version1+"/" + " -v1 " + version1 + " -d2 " + accuracy_data_dir + version2+"/" + " -v2 " + version2)


print('#####################')
print('##  compare precision')
print('#####################')

precision_data_dir = data_dir + "precision/"
precision_out_dir = out_dir + "precision_comparsion/"

if not os.path.isdir(precision_out_dir):
    os.mkdir(precision_out_dir)

os.system("python3 BIP_validation_pipelines_precision/compare_precision_values.py -o " + precision_out_dir + " -d1 " + precision_data_dir + version1+"/" + " -v1 " + version1 + " -d2 " + precision_data_dir + version2+"/" + " -v2 " + version2)
#print("python3 BIP_validation_pipelines_precision/compare_precision_values.py -o " + precision_out_dir + " -d1 " + precision_data_dir + version1+"/" + " -v1 " + version1 + " -d2 " + precision_data_dir + version2+"/" + " -v2 " + version2)


print('#######################')
print('##  compare sensitivity')
print('#######################')

sensitivity_data_dir = data_dir + "sensitivity/"
sensitivity_out_dir = out_dir + "sensitivity_comparsion/"

if not os.path.isdir(sensitivity_out_dir):
    os.mkdir(sensitivity_out_dir)

os.system("python3 BIP_validation_pipelines_LoD/compare_LoD_values.py -i1 " + bip_out_dir + " -i2 " + bip_out_dir + " -o " + sensitivity_out_dir + " -d1 " + sensitivity_data_dir + version1+"/" + " -v1 " + version1 + " -d2 " + sensitivity_data_dir + version2+"/" + " -v2 " + version2)
#print("python3 BIP_validation_pipelines_LoD/compare_LoD_values.py -i1 " + bip_out_dir + " -i2 " + bip_out_dir + " -o " + sensitivity_out_dir + " -d1 " + sensitivity_data_dir + version1+"/" + " -v1 " + version1 + " -d2 " + sensitivity_data_dir + version2+"/" + " -v2 " + version2)


print('#######################')
print('##  compare specificity')
print('#######################')

specificity_data_dir = data_dir + "specificity/"
specificity_out_dir = out_dir + "specificity_comparsion/"

if not os.path.isdir(specificity_out_dir):
    os.mkdir(specificity_out_dir)

os.system("python3 BIP_validation_pipelines_LoB/compare_LoB_values.py -i1 " + bip_out_dir + " -i2 " + bip_out_dir + " -o " + specificity_out_dir + " -d1 " + specificity_data_dir + version1+"/" + " -v1 " + version1 + " -d2 " + specificity_data_dir + version2+"/" + " -v2 " + version2)
#print("python3 BIP_validation_pipelines_LoB/compare_LoB_values.py -i1 " + bip_out_dir + " -i2 " + bip_out_dir + " -o " + specificity_out_dir + " -d1 " + specificity_data_dir + version1+"/" + " -v1 " + version1 + " -d2 " + specificity_data_dir + version2+"/" + " -v2 " + version2)
