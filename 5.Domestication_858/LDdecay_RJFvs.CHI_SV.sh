#!/bin/bash
PopLDdecay -InVCF ./sv.vcf -SubPop ../RJF.txt -MaxDist 2000 -OutStat RJF_LDdecay
PopLDdecay -InVCF ./sv.vcf -SubPop ../Chicken.txt -MaxDist 2000 -OutStat CHI_LDdecay

echo -e "RJF_LDdecay.stat.gz\tRJF" > Pop.list
echo -e "CHI_LDdecay.stat.gz\tCHI" >> Pop.list

perl ~/software/bin/Plot_MultiPop.pl -inList Pop.list -bin1 5000 -bin2 10000 -output RJF_CHI_LDdecay
