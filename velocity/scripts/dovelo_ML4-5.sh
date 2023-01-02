unset Reads
declare -a Reads
for file in $1/*L*R*
do
        Reads+=($file)
done

kb count --h5ad -i ~/tools/Kallisto/Velo_97/cDNA_introns.idx -g ~/tools/Kallisto/Velo_97/tr2g.tsv -x 10xv3 -o $2 \
-c1 ~/tools/Kallisto/Velo_97/cDNA_tx_to_capture.txt -c2 ~/tools/Kallisto/Velo_97/introns_tx_to_capture.txt --lamanno --filter bustools -t 4  ${Reads[@]}