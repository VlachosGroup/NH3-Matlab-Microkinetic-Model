#$ -N amm_main.m
#$ -m eas
#$ -M wittreich@udel.edu
#$ -pe threads 5
#$ -l m_mem_free=1G

vpkg_require matlab/2016b
matlab -nodisplay -r 'amm_main'
