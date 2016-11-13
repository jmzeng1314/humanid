##http://www.genome.jp/kegg-bin/get_htext?hsa00001+3101
## You need to open this link through a browser and then Download the file by click Download htext !
##download files :hsa00001.keg ,about 2.7M (2016/02/16)
perl -lane '{if(/^C/){$t=$F[1]};if(/^D/){print "$t\t$F[1]"}}' hsa00001.keg >kegg2geneID.txt
perl -lane '{if(/^A<b>(.*?)<\/b>/){$a=$1};if(/^B.*?<b>(.*?)<\/b>/){$b=$1};if(/^C\s+\d+\s+(.*?)\s+\[/){print "$a\t$b\t$F[1]\t$1"};}' hsa00001.keg >kegg_hierarchical.txt

