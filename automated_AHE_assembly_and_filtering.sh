##Python scripts require BioPython and some only work with Python 2.7
##All scripts must be placed in your path.
##Scripts: singleline.pl, split.py, removelist.py, s_hit_checker.py, ortholog_filter.py, split.py, extract_probe_region.py, FASconCAT-G_v1.04.pl, counting_monster.py
##A blast database for your reference genome (from the probe design process) needs to be put together. If you have you genome assembly the make blast database command is <makeblastdb -in Biomphalaria.fasta -dbtype nucl> Change Biomphalaria to your genome assembly


#####Execute script in folder containing IBA.py output. I typically run this a folder with only the IBA.py output.


#####VARIABLES THAT YOU NEED TO SPECIFY###########################################################################################################################
REFERENCE="ELithasia_geniculata"  ##Change name to whatever your reference taxon is. The name must be a part of the sequence header for you reference taxon in you reference fasta files. It should be the whole taxon name (e.g., ELithasia_geniculata)
GENOMEFASTADB="/run/media/labgroup/storage/BiomphalariaGenome/biomphalaria_genome.fa"  ##Path to your genome blast database
CORES=31  ##Specify the maximum number of threads/cores to use
PATH_TO_REF_FOLDER="/run/media/labgroup/storage/FWS_123701_REF"  ##This should have you reference alignments that were used to design your probes. Probably 3-6 species.
##################################################################################################################################################################

###This script starts after IBA.py has been run. Please see script processAHE.sh



###Backup output of IBA, which is the longest step so you don't want to lose the original output.
mkdir backup
cp *.fasta backup/
cp *.table backup/


#Concatenate all assembled loci for all taxa into one file
cat *finalseqs.fasta >ALL_FULL_LOCI.fa

##Make fasta file not interleaved
singleline.pl ALL_FULL_LOCI.fa >sl_ALL_FULL_LOCI.fa

##Split sequences by locus. IBA outputs by taxa, but we all taxa for each locus.
split.py sl_ALL_FULL_LOCI.fa _  ## "_" is delimiter...it could be different for your files. For ours it is >LOCUS_TAXON__comp0_seq0. The first underscore is the important part here.

#the resulting files are named by loci  L1.fa-L###.fa in my case. For each locus we used unix command line with a loop to add the new assembled data to an alignment of the reference taxa using mafft and these resulting files will have a prefix of refaa_


##Alignment with mafft using commands that the default executable of mafft didn't seem to have (or I couldn't get it to work). Change path to mafft.bat and the folder with your reference sequences.
for X in L*.fa; do /home/labgroup/programs/mafft-linux64/mafft.bat --auto --maxiterate 1000 --thread $CORES --adjustdirectionaccurately --addlong $X $PATH_TO_REF_FOLDER/$X\s> aligned_$X; done
#TO DO: automate changing path to reference folder.

####CHECKPOINT 1
echo "CHECKPOINT 1"
mkdir checkpoint1
cp aligned_* checkpoint1/


echo "MAFFT is finished. A few steps are proceeding before blast steps"
##Make all alignments single lined
for X in aligned_L* ; do singleline.pl $X > 1l_$X; done

##Create Text File with a list of your aligned sequences to be used in next step.
ls 1l_aligned_L* >list.txt

###create large file with all sequences and taxa, including refernces, for use later
cat 1l_aligned_L*.fa >ALL_FULL_LOCI_REF.fa
sed -i 's/-//g' ALL_FULL_LOCI_REF.fa
sed -i 's/>_R_/>/' ALL_FULL_LOCI_REF.fa

##Extract Probe, head, and tail regions
extract_probe_region.py list.txt $REFERENCE outdir ###Biomphalaria is the reference taxon used to design probes...it should be the same as what you used in IBA.py ##outdir is the director to output the file. Can be renamed
cd outdir
sed -i 's/>_R_/>/' 1l_aligned_L* ##MAFFT puts a _R_ for sequences it reversed during alignment. This should be removed
cat *.trimmed > 3by3.in



##Blast all sequences to check for single hits and remove sequences with greater than 1 hit.
sed -i 's/-//g' 3by3.in  ##First remove gaps as blast does not like some (or all) of the alignments

echo "3by3 blast is running"
tblastx -query 3by3.in -db $GENOMEFASTADB -out 3by3.out -outfmt 6 -max_target_seqs 3 -max_hsps 3 -num_threads $CORES ##Can also do tblastx if sequences are divergent from reference.
echo "3by3 blast is finished"

s_hit_checker.py 3by3.out 0.90
removelist.py 3by3.in 3by3_del_list0.90.txt 1by1.in


#Filter for sequences that passed single hit criteria and found best hit with blast and filtered for orthologs homology
echo "1by1 blast and filtering is starting"
tblastx -query 1by1.in -db $GENOMEFASTADB -out 1by1.out -outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads $CORES  ##Can also do tblastx if sequences are divergent from reference.
ortholog_filter.py 1by1.out $REFERENCE
getlist.py 1by1.in 1by1_keep_list ORTHO_PASS.fa

echo "1by1 blast and filtering for non-orthologs is finished"

###CHECKPOINT 2
echo "CHECKPOINT 2"
mkdir checkpoint2/
cp ORTHO_PASS.fa checkpoint2/


##Make single line fasta (not interleaved) and split sequences that passed orthology into single locus files.
singleline.pl ORTHO_PASS.fa >sl_ORTHO_PASS.fa

sed -i "s/_/|/g" sl_ORTHO_PASS.fa
sed -i "s/|seq/_seq/g" sl_ORTHO_PASS.fa
split.py sl_ORTHO_PASS.fa \|
for X in L*.fa; do /home/labgroup/programs/mafft-linux64/mafft.bat --thread $CORES $X > $X\s; done

##collapse heterogeneity into a consensus for each taxon and locus
FASconCAT-G_v1.04.pl -c -c -c -o -s
cat FcC_L* > NOISO_PROBE.fa  ##No isoforms, just consensus with IUPAC ambuguity codes for DNA
sed -i "s/|/_/g" NOISO_PROBE.fa
sed -i "s/_consensus//" NOISO_PROBE.fa
remove_duplicates.py NOISO_PROBE.fa


###CHECKPOINT 3
echo "Checkpoint 3"
mkdir checkpoint3
cp NOISO_PROBE.fa checkpoint3/

###use usearch to get sequences using a feature that allows a partial match since we have removed part of the original sequence name
usearch -fastx_getseqs ../ALL_FULL_LOCI_REF.fa -labels NOISO_PROBE_keep.list -label_substr_match -fastaout FINAL_FULL_CLEAN_REF.fa
sed -i 's/>_R_/>/' FINAL_FULL_CLEAN_REF.fa  ##Get rid of _R_ added by MAFFT when making alignments for ALL_FULL_LOCI_REF.fa
singleline.pl FINAL_FULL_CLEAN_REF.fa > 1l_FINAL_FULL_CLEAN_REF.fa  #make single line, not interleaved.
sed -i "s/_/|/g" 1l_FINAL_FULL_CLEAN_REF.fa ##Clean up headers
sed -i "s/|seq/_seq/g" 1l_FINAL_FULL_CLEAN_REF.fa ##Clean up some more

mkdir gapA
mv 1l_FINAL_FULL_CLEAN_REF.fa gapA/
cd gapA
split.py 1l_FINAL_FULL_CLEAN_REF.fa \|  ##Split into individual loci




###CHECKPOINT 4
echo "CHECKPOINT 4"
mkdir checkpoint4
cp L*.fa checkpoint4/

###Align again
for X in L*.fa; do /home/labgroup/programs/mafft-linux64/mafft.bat --thread $CORES --adjustdirectionaccurately --allowshift --unalignlevel 0.8 --leavegappyregion --maxiterate 1000 --globalpair $X > gapA_$X\s; done
sed -i 's/>_R_/>/' gapA_L*.fas  ##Get rid of _R_ to start some fasta headers
mkdir nextStepAlignments
cp gapA_L*.fas nextStepAlignments/
cd nextStepAlignments/
FASconCAT-G_v1.04.pl -c -c -c -o -s  ##Collapse heterogeneity into a consensus using ambigious IUPAC codes. May want to explore for coalescent based analyses.


cat FcC_gapA_L* > ALL_single_for_CM.fas ##file that can be used with counting_monster.pyt
singleline.pl  ALL_single_for_CM.fas > sl_ALL_single_for_CM.fas
sed -i 's/_consensus//' sl_ALL_single_for_CM.fas
counting_monster.py sl_ALL_single_for_CM.fas \|  ##counting monster makes a table (tableout.txt) that has which taxa are present in which loci with total numbers at end of columns and rows. Useful and quick data evaluation.
sed -i "s/|/_/g" FcC_gapA_L*
sed -i -r "s/_+comp.\+//" FcC_gapA_L*
sed -i "s/^>L[0-9]\+_/>/" FcC_gapA_L*
sed -i 's/>_/>/' FcC_gapA_L*
sed -i 's/_R/__comp0/' FcC_gapA_L*
for X in FcC_*; do singleline.pl $X > 1l_$X; done
ls 1l_FcC* > list.txt
extract_probe_region.py list.txt $REFERENCE outdir ##This is where probe region is seperated from sequence fragments before and after the probe region. The head and tail regions are captured for many taxa during the target capture process.
cd outdir
rename .fas.trimmed .fas *.fas.trimmed
rename FcC_gapA_L L FcC_gapA_L*

###Move Final probe region alignments to folder in starting directory
mkdir ../../../../probeRegionAlignments_FINAL
mv *.fas ../../../../probeRegionAlignments_FINAL/



###Alignments will be in probeRegionAlignments folder. If you want sequence regions before and after probe regions, you will need to make some modifications to concatenate them or not split them apart.

