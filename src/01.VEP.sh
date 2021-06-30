SCRATCH=~/scratch/09.COVID19/05.BQC/08.OAS1/
ANNO_DIR=/home/richards/tomoko.nakanishi/my_project/13.UKB_WES/scratch/200KVEP

conda activate vep

#ENSG00000174840 PDE12 pLI/LOEUF = 0.04/0.6
chr=3
filter_vep -i ${ANNO_DIR}/WES_200K.${chr}.var.lof.cadd.gnomad.txt -o stdout -filter "CADD_PHRED > 20 and Gene is ENSG00000174840" | grep -v "#" | cut -f1 | sort | uniq | awk '{print "chr"$0}' > ${SCRATCH}/ENSG00000174840.CADD.list &
filter_vep -i ${ANNO_DIR}/WES_200K.${chr}.var.lof.cadd.gnomad.txt -o stdout -filter "LoF is HC and Gene is ENSG00000174840" | grep -v "#" | cut -f1 | sort | uniq | awk '{print "chr"$0}' > ${SCRATCH}/ENSG00000174840.LoF.list &

#ENSG00000118507 AKAP7 pLI/LOUEF = 0/1.22
chr=6
filter_vep -i ${ANNO_DIR}/WES_200K.${chr}.var.lof.cadd.gnomad.txt -o stdout -filter "CADD_PHRED > 20 and Gene is ENSG00000118507" | grep -v "#" | cut -f1 | sort | uniq | awk '{print "chr"$0}' > ${SCRATCH}/ENSG00000118507.CADD.list &
filter_vep -i ${ANNO_DIR}/WES_200K.${chr}.var.lof.cadd.gnomad.txt -o stdout -filter "LoF is HC and Gene is ENSG00000118507" | grep -v "#" | cut -f1 | sort | uniq | awk '{print "chr"$0}' > ${SCRATCH}/ENSG00000118507.LoF.list &

