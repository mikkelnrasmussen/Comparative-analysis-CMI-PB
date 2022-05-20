# Created by Joaquin Reyna, source: https://github.com/joreynajr/cmi-pb-multiomics

# set bash strict mode 
set -euo pipefail
IFS=$'\n\t'

outdir="raw/"
mkdir -p $outdir

#########################################################################################
# Download 2020 data ####################################################################
#########################################################################################
#2020LD_subject.csv
#2020LD_specimen.csv
#2020LD_live_cell_percentages.csv
#2020LD_olink_prot_exp.csv
#2020LD_rnaseq.csv

wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2020LD_subject.csv
wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2020LD_specimen.csv
wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2020LD_live_cell_percentages.csv
wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2020LD_olink_prot_exp.csv
wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2020LD_rnaseq.csv

#########################################################################################
# Download 2021 data ####################################################################
#########################################################################################
#2021BD_subject.csv
#2021BD_specimen.csv
#2021BD_live_cell_percentages.csv
#2021BD_olink_prot_exp.csv
#2021BD_rnaseq.csv

wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2021BD_subject.csv
wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2021BD_specimen.csv
wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2021BD_live_cell_percentages.csv
wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2021BD_olink_prot_exp.csv
wget -P $outdir https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/2021_cmipb_challenge/2021BD_rnaseq.csv

