# $1: path to fastq files
# $2: path to output of kallisto
# $3: path to barcodes

unset Reads
declare -a Reads
for file in $1/*R*
do
	Reads+=($file)
done

kallisto bus -i ~/tools/Kallisto/Velo/cDNA_introns.idx -o $2 -x 10xv2 -t6 ${Reads[@]}

bustools correct -w ~/tools/Kallisto/whitelist_v2.txt -p $2/output.bus | bustools sort -o $2/output.correct.sort.bus -t4 -
bustools capture -s -x -o $2/spliced.bus -c ~/tools/Kallisto/Velo/introns_tx_to_capture.txt -e $2/matrix.ec -t $2/transcripts.txt $2/output.correct.sort.bus
bustools capture -s -x -o $2/unspliced.bus -c ~/tools/Kallisto/Velo/cDNA_tx_to_capture.txt -e $2/matrix.ec -t $2/transcripts.txt $2/output.correct.sort.bus

bustools count -o $2/unspliced -g ~/tools/Kallisto/Velo/tr2g.tsv -e $2/matrix.ec -t $2/transcripts.txt --genecounts $2/unspliced.bus
bustools count -o $2/spliced -g ~/tools/Kallisto/Velo/tr2g.tsv -e $2/matrix.ec -t $2/transcripts.txt --genecounts $2/spliced.bus
unset Reads
echo done
python3 scripts/KallistoLoom.py $2 $3