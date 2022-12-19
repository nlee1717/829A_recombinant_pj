variant='SA'
ref_seq=/main/dir/ref_seq.fasta
multiR=/main/dir/SA_find_multiallelic.R
cd /main/dir/$variant
# there should be separate folders for each variants in the /main/dir 
echo $'\n##### Running pipeline for '$variant$'#####\n\n'

# head -10 $variant'_data.fasta'

# gets GISAID ids for reproduction purpose
grep ">" $variant'_data.fasta' | awk -F"|" '{ print $2 }' > $variant'_id.txt';
# head $variant'_id.txt'

echo $'\n\n========== Alignment =========='
echo $'\n<< mafft alignment >>'
mafft --auto --thread -1 --keeplength --quiet --mapout --preservecase --addfragments $variant'_data.fasta' $ref_seq > $variant'_alignment.fasta'

echo $'\n<< replace symbols with N: '$variant'_alignment_cleaned.fasta>>'
sed '/^>/! s/[^actgACTG]/N/g' $variant'_alignment.fasta' > $variant'_alignment_cleaned.fasta'

echo $'\n<< removing some seqs: '$variant'_alignment_good.fasta >>'
grep "Stretches" $variant'_sequencing_info.tsv' | awk '{print $1 " " $2 "|" $3 "|" $4}' > Ns_id.txt;
sed -i.bak "s/\ /_/g" Ns_id.txt;
sed -i.bak "s/\ /_/g" $variant'_alignment_cleaned.fasta';
seqkit grep -n -v -f Ns_id.txt $variant'_alignment_cleaned.fasta' > $variant'_alignment_good.fasta'

echo $'\n<< check size with seqkit >>'
seqkit stats $variant'_alignment_good.fasta'

echo $'\n\n========== Split by Strain =========='
echo $'\n<< extract specific strain: '$variant'_newstrain_id.txt >>'
awk '{print $1 " " $2 "|" $3 "|" $4}' $variant'_newstrain_sequencing_info.tsv' | grep "hCoV" > $variant'_newstrain_id.txt';
sed -i.bak "s/\ /|/g" $variant'_newstrain_id.txt'
sed -i.bak "s/|Asia//" $variant'_newstrain_id.txt'

echo $'\n !!! CHECK IF PATTERNS MATCH !!!\n'
echo '> from newstrain_id:'
head -5 $variant'_newstrain_id.txt'
echo $'\n> from aligned fasta:'
head -1 $variant'_alignment_good.fasta'
echo $'\n'

while true; do
    read -p ">>> Press n to stop. Any other key to continue " ans
    case $ans in
        [n]* ) exit;;
        * ) break;;
    esac
done

seqkit grep -n -f $variant'_newstrain_id.txt' $variant'_alignment_good.fasta' > $variant'_newstrain.fasta';
seqkit grep -n -v -f $variant'_newstrain_id.txt' $variant'_alignment_good.fasta' > $variant'_oldstrain.fasta'

echo $'\n<< check exact dups: '$variant'_XX_nondup.fasta >>'
seqkit rmdup -D$variant'_olddeleted.txt' -s $variant'_oldstrain.fasta' > $variant'_oldstrain_nondup.fasta';
seqkit rmdup -D$variant'_newdeleted.txt' -s $variant'_newstrain.fasta' > $variant'_newstrain_nondup.fasta';
rm $variant'_newstrain.fasta';
rm $variant'_oldstrain.fasta'

echo $'\n<< check size with seqkit >>'
seqkit stats $variant'_oldstrain_nondup.fasta';
seqkit stats $variant'_newstrain_nondup.fasta'

echo $'\n<< remove ref seq from : '$variant'_oldstrain_nondup.fasta >>'
grep -v 'EPI_ISL_402125' $variant'_oldstrain_nondup.fasta' > temp.fasta
sed -n '/>/,$p' temp.fasta > $variant'_oldstrain_nondup.fasta'
rm temp.fasta
head -1 $variant'_oldstrain_nondup.fasta'

echo $'\n\n========== Sample Generation =========='
echo $'\n<< build sample: '$variant'_sample.fasta>>'
cat $ref_seq > $variant'_sample.fasta';
seqkit sample -s 43984291 -p 0.6 $variant'_oldstrain_nondup.fasta' | seqkit shuffle -s 92834717 | seqkit head -n 25 >> $variant'_sample.fasta';
seqkit sample -s 23849817 -p 0.6 $variant'_newstrain_nondup.fasta' | seqkit shuffle -s 34876261 | seqkit head -n 25 >> $variant'_sample.fasta';

echo $'\n<< get accession numbers: '$variant'_sample_id.txt >>'
grep ">" $variant'_sample.fasta' | awk -F"|" '{ print $2 }' > $variant'_sample_id.txt';
cat $variant'_sample_id.txt'

echo $'\n\n========== Masking problematic sites =========='
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $variant'_sample.fasta' > $variant'_sample_temp.fasta';
sed 's/./& /g' < $variant'_sample_temp.fasta' > $variant'_sample_test.fasta';
rm $variant'_sample_temp.fasta';
Rscript $multiR $variant'_sample_test.fasta' $variant'_sample_masked.fasta' $variant'_sample_positions.txt'

cat $variant'_sample_positions.txt'

echo $'\n\n========== Labelling =========='
seqkit head -n 1 $variant'_sample_masked.fasta' > $variant'_sample_masked_id.fasta';
seqkit range -r 2:26 $variant'_sample_masked.fasta' > f1.f;
seqkit range -r 27:51 $variant'_sample_masked.fasta' > f2.f;
awk '/^>/{print ">O" ++i; next}{print}' < f1.f >> $variant'_sample_masked_id.fasta';
awk '/^>/{print ">N" ++i; next}{print}' < f2.f >> $variant'_sample_masked_id.fasta';
rm f1.f;
rm f2.f

echo $'\n\n========== KwARG =========='
/mnt/c/Users/nklee/Desktop/binaries_linux/kwarg -T50,30 -Q500 -S-1,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,1 -M-1,1.91,1.81,1.71,1.61,1.51,1.41,1.31,1.21,1.11,1.01,0.91,0.81,0.71,0.61,0.51,0.41,0.31,0.21,0.11,0.02,1.1 -R1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1 -C2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,-1 -k -n -f < $variant'_sample_masked_id.fasta' > $variant'_kwarg_out.txt'
# /mnt/c/Users/nklee/Desktop/binaries_linux/kwarg -Q100 -S-1,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,1 -M-1,1.91,1.81,1.71,1.61,1.51,1.41,1.31,1.21,1.11,1.01,0.91,0.81,0.71,0.61,0.51,0.41,0.31,0.21,0.11,0.02,1.1 -R1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1 -C2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,-1 -k -n -f < delta_sample_masked_id.fasta > delta_kwarg_out.txt
# /mnt/c/Users/nklee/Desktop/binaries_linux/kwarg -T50,30 -Q500 -S-1,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,1 -M-1,1.91,1.81,1.71,1.61,1.51,1.41,1.31,1.21,1.11,1.01,0.91,0.81,0.71,0.61,0.51,0.41,0.31,0.21,0.11,0.02,1.1 -R1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1 -C2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,-1 -k -n -f -Z -dtest_arg.dot -e < $variant'_sample_masked_id.fasta'
# # add seed number after -Z
