#!/bin/bash
GTF_FILE="reference_ensembl_org_pub_grch37_release-84_fasta_homo_sapiens/Homo_sapiens.GRCh37.82.chr.gtf"
BAM_DIR="bams"
N_BAMS=
T_BAMS=
for i in {1..5}
do
   if [ ${i} -ne 1 ]
   then
      N_BAMS="${N_BAMS},"
      T_BAMS="${T_BAMS},"
   fi
   N_BAMS="${N_BAMS}${BAM_DIR}/${i}N/tophat/accepted_hits.bam"
   T_BAMS="${T_BAMS}${BAM_DIR}/${i}T/tophat/accepted_hits.bam"
done

CMD="cuffdiff -p 1 ${GTF_FILE} ${N_BAMS} ${T_BAMS}"
echo $CMD
