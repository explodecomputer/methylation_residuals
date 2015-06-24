#!/bin/bash

kidsgeno="${HOME}/repo/mQTL-partitioning/grm/alspac_hm3_kids_un"
mothersgeno="${HOME}/repo/mQTL-partitioning/grm/alspac_hm3_mothers"

kidsids="${HOME}/repo/methylation_residuals/data/15up.id"
mothersids="${HOME}/repo/methylation_residuals/data/FOM.id"

kidsout="${HOME}/repo/methylation_residuals/grms/kids"
mothersout="${HOME}/repo/methylation_residuals/grms/mothers"

kidsmgrm="${kidsout}.mgrm"
mothersmgrm="${mothersout}.mgrm"

rm -f ${kidsmgrm}
rm -f ${mothersmgrm}


for i in {1..22}
do
	plink1.90 --bfile ${kidsgeno} --make-grm-bin --keep ${kidsids} --maf 0.01 --chr ${i} --out ${kidsout}${i}
	echo ${kidsout}${i} >> ${kidsmgrm}
done

for i in {1..22}
do
	plink1.90 --bfile ${mothersgeno} --make-grm-bin --keep ${mothersids} --maf 0.01 --chr ${i} --out ${mothersout}${i}
	echo ${mothersout}${i} >> ${mothersmgrm}
done


gcta64 --mgrm ${kidsmgrm} --make-grm --out ${kidsout}
gcta64 --mgrm ${mothersmgrm} --make-grm --out ${mothersout}

gcta64 --grm ${kidsout} --grm-cutoff 0.05 --make-grm --out ${kidsout}
gcta64 --grm ${mothersout} --grm-cutoff 0.05 --make-grm --out ${mothersout}


for i in {1..22}
do
	gcta64 --grm ${kidsout}${i} --keep ${kidsout}.grm.id --make-grm --out ${kidsout}${i}
	echo ${kidsout}${i} > ${kidsout}${i}.mgrm
	echo ${kidsout} >> ${kidsout}${i}.mgrm
done

for i in {1..22}
do
	gcta64 --grm ${mothersout}${i} --keep ${mothersout}.grm.id --make-grm --out ${mothersout}${i}
	echo ${mothersout}${i} > ${mothersout}${i}.mgrm
	echo ${mothersout} >> ${mothersout}${i}.mgrm
done

