#!/bin/bash

#PBS -V
#PBS -m e -M ntz@stsci.edu
#SBATCH --job-name=APLC_quart_hex1X025c1_N0250_FPM400M050_LSann30D85X025c1Pad05_Img100C_35WA100CA100_BW10Nlam03fpres2_linbarpre1
#SBATCH -o ./nccs_survey_test/logs/APLC_quart_hex1X025c1_N0250_FPM400M050_LSann30D85X025c1Pad05_Img100C_35WA100CA100_BW10Nlam03fpres2_linbarpre1.log
#SBATCH --account=s1649

#SBATCH --constraint=hasw
#SBATCH --ntasks=1 --nodes=1

#SBATCH --qos=allnccs
#SBATCH --time=12:00:00

. /usr/share/modules/init/bash
module purge
module load comp/intel-10.1.017
ulimit -s unlimited

#Optional: monitor the memory usage...
mkdir -p ${NOBACKUP}/policeme
/usr/local/other/policeme/policeme.exe -d ${NOBACKUP}/policeme

ampl ./nccs_survey_test/amplsrc/APLC_quart_hex1X025c1_N0250_FPM400M050_LSann30D85X025c1Pad05_Img100C_35WA100CA100_BW10Nlam03fpres2_linbarpre1.mod

exit 0