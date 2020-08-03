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

    2. Run LoD scripts over the different samples and different bip versions

    3. create an LoD summary table for the variants of interest in the study

    4. create a reportable table comparing the LoD values between the two versions of bip
    """

parser = argparse.ArgumentParser(description="BIP LoD comparison runner",epilog=epilog, formatter_class=CustomFormatter)
parser.add_argument("-i1", "--in_dir1",
                    action="store",
                    dest="in_dir1",
                    required=True,
                    help=("the directory where the BIP(any version) outputs are stored"))
parser.add_argument("-i2", "--in_dir2",
                    action="store",
                    dest="in_dir2",
                    required=True,
                    help=("the directory where the BIP(any version) outputs are stored"))
parser.add_argument("-o", "--output_dir",
                    action="store",
                    dest="output_dir",
                    required=True,
                    help=("directory to store the tables"))
parser.add_argument("-d1", "--_data_tables_dir1",
                    action="store",
                    dest="_data_tables_dir1",
                    required=True,
                    help=("directory of the sample and variant tables for the study"))
parser.add_argument("-v1", "--version1",
                    action="store",
                    dest="version1",
                    required=True,
                    help=("specify the version of the BIP used"))
parser.add_argument("-d2", "--_data_tables_dir2",
                    action="store",
                    dest="_data_tables_dir2",
                    required=True,
                    help=("directory of the sample and variant tables for the study"))
parser.add_argument("-v2", "--version2",
                    action="store",
                    dest="version2",
                    required=True,
                    help=("specify the version of the BIP used"))


args = parser.parse_args()

if not os.path.isdir(args.in_dir1) and not os.path.isdir(args.in_dir2):
    raise RuntimeError("Nonexistant or invalid input directory '{}'".format(args.in_dir))

abs_dir = subprocess.Popen(["pwd"], stdout=subprocess.PIPE, shell=True).communicate()[0].decode('utf-8').replace('\n','/')

docker_image="docker-scratch.artifactory01.ghdna.io/rgh-mosesdev:1.0.0"
print(abs_dir)
r_script_dir = abs_dir
r_script="D000630_lod_report.MAIN.200128_custom.R"
bip_out_dir=args.in_dir1
data_tables_dir=args._data_tables_dir1
out_dir=args.output_dir
version1=args.version1

#file_dir = os.path.dirname(os.path.realpath(__file__))
#print(abs_dir)
#print(file_dir)

out_dir1 = out_dir+'/out_1_bip_'+version1+'/'
if not os.path.isdir(out_dir1):
    os.makedirs(out_dir1)
#os.system("docker exec --user 36154:21508 ade6b00f8f33 Rscript " + abs_dir + r_script + " -i " + bip_out_dir + " -d " + data_tables_dir + " -o " + out_dir1 + " -v " + version1 + " -a " + abs_dir)

os.system("docker run --rm -v /ghds/:/ghds/ --user $(id -u):$(id -g) " + docker_image + " Rscript " + r_script_dir + r_script + " -i " +  bip_out_dir + " -d " + data_tables_dir + " -o " +out_dir1 + " -v " + version1 + " -a " + abs_dir)

os.system("python3 makeLoDSummaryTable.py " + out_dir1 + " " + version1)

bip_out_dir=args.in_dir2
data_tables_dir=args._data_tables_dir2
version2=args.version2

out_dir2 = out_dir+'/out_2_bip_'+version1+'/'
if not os.path.isdir(out_dir2):
    os.makedirs(out_dir2)
#os.system("docker exec --user 36154:21508 ade6b00f8f33 Rscript " + abs_dir + r_script + " -i " + bip_out_dir + " -d " + data_tables_dir + " -o " + out_dir2 + " -v " + version2 + " -a " + abs_dir)

os.system("docker run --rm -v /ghds/:/ghds/ --user $(id -u):$(id -g) " + docker_image + " Rscript " + r_script_dir + r_script + " -i " +  bip_out_dir + " -d " + data_tables_dir + " -o " +out_dir2 + " -v " + version2 + " -a " + abs_dir)

os.system("python3 makeLoDSummaryTable.py " + out_dir2 + " " + version2)

os.system("python3 combineLoDTables.py " + out_dir1 + " " + out_dir2 + " " + version1 + " " + version2 + " " + out_dir)
