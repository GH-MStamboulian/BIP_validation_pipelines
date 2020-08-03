#!/bin/bash

docker_image=docker-scratch.artifactory01.ghdna.io/rgh-mosesdev:1.0.0

#abs_dir = subprocess.Popen(["pwd"], stdout=subprocess.PIPE, shell=True).communicate()[0].decode('utf-8').replace('\n','/')
abs_dir=$(pwd)"/"

r_script_dir=$abs_dir
r_script="D000630_lod_report.MAIN.200128_custom.R"

bip_out_dir="/ghds/ivd/flowcentral/"
data_tables_dir="/ghds/groups/bioinformatics/02_DEVELOPMENT/moses/D000630_data_bip_3_5_3/"
out_dir="/ghds/groups/bioinformatics/02_DEVELOPMENT/moses/D000630_bip_3_5_3_docker_results/"
version="3.5.3"

echo ${r_script_dir}${r_script}

docker run --rm -v /ghds/:/ghds/ --user $(id -u):$(id -g) ${docker_image} Rscript ${r_script_dir}${r_script} -i ${bip_out_dir} -d ${data_tables_dir} -o ${out_dir} -v ${version} -a ${abs_dir}
