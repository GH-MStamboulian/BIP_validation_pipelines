#!/bin/bash

module load sge/2011.11p1
module load python/3.6.5

flowCell_f="flowcells.txt"
sampleSheet_dir="/ghds/ivd/flowcentral/"
out_dir=""
#n_cores=5
docker_img="docker.artifactory01.ghdna.io/bip:3.5.3"
panel="GH2.11"
raw_dir="/ghds/ivd/raw/"
suffix=".3.5.3_rerun"

parent_dir=$(pwd)

while IFS= read -r flowCellID
do
    out_folder=${flowCellID%"$suffix"}
    if [ ! -d $out_dir$out_folder ]; then
	mkdir $out_dir$out_folder
        cp -r $sampleSheet_dir$flowCellID"/SampleSheet_new.csv" $out_dir$out_folder"/SampleSheet.csv"
	python3 /ghds/shared/apps/bip_sge_launcher/mono-bip.py -d $docker_img --outdir $parent_dir"/"$out_dir$out_folder"/" --sample_sheet $parent_dir"/"$out_dir$out_folder"/SampleSheet.csv" --parameter_set $panel --seqdir $raw_dir$out_folder"/" --launch
    fi
	#echo python3 /ghds/shared/apps/bip_sge_launcher/mono-bip.py -t $n_cores -d $docker_img --outdir $parent_dir"/"$out_dir$out_folder"/" --sample_sheet $parent_dir"/"$out_dir$out_folder"/SampleSheet.csv" --parameter_set $panel --seqdir $raw_dir$out_folder"/" --launch
    	
done < "$flowCell_f"

#python3 /ghds/shared/apps/bip_sge_launcher/mono-bip.py -t 75 -d docker.artifactory01.ghdna.io/bip:3.5.3 --outdir /ghds/groups/bioinformatics/02_DEVELOPMENT/moses/bip_3_5_3/181009_NB551348_0023_AHKG22BGX7/ --sample_sheet /ghds/groups/bioinformatics/02_DEVELOPMENT/moses/bip_3_5_3/181009_NB551348_0023_AHKG22BGX7/SampleSheet.csv --parameter_set GH2.11 --seqdir /ghds/ivd/raw/181009_NB551348_0023_AHKG22BGX7/ --launch

#python3 /ghds/shared/apps/bip_sge_launcher/mono-bip.py -t 5 -d docker.artifactory01.ghdna.io/bip:3.5.3 --outdir 181221_NB551348_0044_AHF5TMBGX7 --sample_sheet 181221_NB551348_0044_AHF5TMBGX7/SampleSheet.csv --parameter_set GH2.11 --seqdir /ghds/ivd/raw/181221_NB551348_0044_AHF5TMBGX7/ --launch
