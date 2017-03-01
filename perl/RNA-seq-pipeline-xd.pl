use strict;
# SGE_Batch -c ' perl /capecchi/pharmacy/morgunlab/Dong/projects/CVID/2016.1.RNA-SEQ.humanCellCulture/step1-workFlow-cutadapt.fq.pl $SGE_TASK_ID' -t 1-13 -P 8 -f 50G -F 50G -r cvidHat2
my $samtoolsBin = "/capecchi/pharmacy/morgunlab/Dong/softwares/RNA-SEQ/samtools/bin/samtools";
my $inputFolder = "/nfs1/Morgun_Lab/ourSequencingData/RNASEQ/2016.1.HumanCellCulture/rawReads";
my $cutadaptTemplate = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC -m 60 --max-n=1 -q 20,20 -o output.fastq input.fastq.gz";

BEGIN{unshift @INC, "/capecchi/pharmacy/morgunlab/Dong/softwares/module/"};
my $softwaresFoldePath = "/capecchi/pharmacy/morgunlab/Dong/softwares";
my $gf = "/capecchi/pharmacy/morgunlab/Dong/database/genomes/human.tophat/Homo_sapiens.GRCh38.79.gtf";
my $bowind = "/capecchi/pharmacy/morgunlab/Dong/database/genomes/human.tophat/Homo_sapiens.GRCh38.dna.toplevel.sm";
my $tophatProgram = "/capecchi/pharmacy/morgunlab/Dong/softwares/RNA-SEQ/tophat/tophat-2.0.13.Linux_x86_64/tophat2";


#unshift @INC, "$softwaresFoldePath/module";
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	#	work flow
	#	1. remove adapters
	#	2. remove low quality sequences(better performed before remove of "N",for some file ,every sequence has an "N" in a particular possition, if we remove "N" first, all sequences will be removed. but these "N"s may be removed in this quality filter step, so after it is removed, we do not have to worry about all sequences to be discarded)
	#	3. remove sequence with any "N"
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 														remove sequence with "N"====record change of file after each step
#																	use prinseq
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

my $outFolder = "$inputFolder.afterCutAdapt";
`mkdir $outFolder` if ( ! -d $outFolder);

my $recordFilteredFile = "$outFolder/recordFilteredSequence_Prinseq.tsv";
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


sub cutadapt{
	my ($inputFile, $OutPutFile) = @_;
	# input 
		# gziped fastq file(string)
		# 3. trimmomatic commond line commond(string) 
	# output
		# trimmomatic result gziped fastq file
#================================
	$cutadaptTemplate =~s/input.fastq.gz/$inputFile/g;
	$cutadaptTemplate =~s/output.fastq/$OutPutFile/g;
	print $cutadaptTemplate,"\n\n"; #die;
	`$cutadaptTemplate` ;
}
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# SGE_Batch -c ' perl /capecchi/pharmacy/morgunlab/Dong/projects/diabetes/RNA-SEQ/2016.1.IRGM1.KO.mice/step1-workFlow-cutadapt.fq.pl $SGE_TASK_ID' -t 1-17 -f 10G -F 10G -r irgmHat

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# do gz files
use sequencePreprocessing::FileHandle;


sub processFastq{
	my @gzFiles = `ls $inputFolder | grep 'fastq' `;
	#my @gzFiles = FileHandle->FindFileWithAppendix($inputFolder,'fastq');
	#print @gzFiles;die;
	#unshift @gzFiles, "/projects1/antibiotic_mice/dong/test/testPrinseq/testTrimmomatic.fastq.gz";	
#	print scalar ( @gzFiles ); die;
	if ($ARGV[0] =~m/[0-9]/){
		my $FileTodo = $ARGV[0]-1;
		@gzFiles = @gzFiles[$FileTodo];
	}
	print @gzFiles,"\nhere\n"; #die;
	foreach my $fileName (@gzFiles)  {
			chomp $fileName;
			print "$fileName\n"; # next;
			# run prinseq, record the number of filtered sequenses
			my $InputFastqFile = "$inputFolder/$fileName";
			my $OutPutFile = "$outFolder/$fileName.cutadapt";
			
			cutadapt($InputFastqFile, $OutPutFile) if (! -e "$OutPutFile");
			my $cutAdaptorFile = $OutPutFile;
			my $tophatOutFolder = "$outFolder.tophat.SGE/$fileName";
			`mkdir $outFolder.tophat.SGE`;
			runTopHat($cutAdaptorFile,$tophatOutFolder) if (! -e "$tophatOutFolder/accepted_hits.bam");
			#`rm $InputFastqFile`;
#			die;	
			my $samFolder = "$outFolder.tophat.SGE.Aligned.sam";
			`mkdir $samFolder` if (! -d $samFolder);
			my $acceptedBam = "$tophatOutFolder/accepted_hits.bam";
			# sort by name, convert to SAM for htseq-count
			print "$samtoolsBin sort -n   $acceptedBam  \n";
			`$samtoolsBin sort -n   $acceptedBam  $acceptedBam.sn ` if (! -e "$samFolder/$fileName");
			`$samtoolsBin  view -o   $samFolder/$fileName   $acceptedBam.sn.bam` if (! -e "$samFolder/$fileName");
			
			my $countFolder = "$outFolder.tophat.SGE.Aligned.sam.htseqCount";
			`mkdir $countFolder` if (! -d $countFolder);
			my $htseqCommand = "/capecchi/pharmacy/morgunlab/Dong/softwares/HTSeq-0.6.1/scripts/htseq-count -s no -a 10 --stranded=no   $samFolder/$fileName  $gf  > $countFolder/$fileName";
			print "\n\n$htseqCommand\n\n\n"; #die;
			`$htseqCommand\n` if (! -e "$countFolder/$fileName");
		} 
}
sub runTopHat {
	my ($cutAdaptorFile, $tophatOutFolder) = @_;
	my $cmd = "$tophatProgram --no-coverage-search -G  $gf  -p 8 -o $tophatOutFolder $bowind $cutAdaptorFile";
	print ($cmd,"\n\n");
	`$cmd`;
	
	
}
processFastq;
#processGZ;
#SGE_Batch -c ' perl workFlow.pl $SGE_TASK_ID ' -t 1-2 -f 4G -F 10G -r copy_StdOuts 														
																
																
																
																
																
																
				
