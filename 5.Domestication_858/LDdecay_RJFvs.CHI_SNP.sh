#!/bin/bash
PopLDdecay -InVCF ./snp.vcf -SubPop ./RJF.txt -MaxDist 500 -OutStat RJF_LDdecay
PopLDdecay -InVCF ./snp.vcf -SubPop ./Chicken.txt -MaxDist 500 -OutStat CHI_LDdecay

echo -e "RJF_LDdecay.stat.gz\tRJF" > Pop.list
echo -e "CHI_LDdecay.stat.gz\tCHI" >> Pop.list

perl ~/software/bin/Plot_MultiPop.pl -inList Pop.list -bin1 10 -bin2 500 -output RJF_CHI_LDdecay
