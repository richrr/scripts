use strict;
use POSIX qw(strftime);
# perl -MGetopt::Long -e 'print( "got it\n" )'
use warnings;
use Getopt::Long qw(GetOptions);
 
my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;

# print the versions of the software used
`echo "$now_string \n\ncutadapt --version" > log.txt` ;
`cutadapt --version >> log.txt` ;

`echo "\nbowtie2 --version" >> log.txt` ;
`bowtie2 --version >> log.txt` ;

`echo "\ntophat -v" >> log.txt` ;
`tophat -v >> log.txt` ;

`echo "\nsamtools --version" >> log.txt` ;
`samtools --version >> log.txt` ;

`echo "\nhtseq-count --help | tail -n 3" >> log.txt` ;
`htseq-count --help | tail -n 3 >> log.txt` ;


# Quality trimming happens before adapter removal http://cutadapt.readthedocs.io/en/latest/guide.html#modifying-reads
# --max-n Discard reads containing more than COUNT N bases. 
# -q trim the 5' and 3' ends with two comma-separated cutoffs (e.g. 15,10)
# -m --minimum-length
# my $cutadaptTemplate = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC -m 60 --max-n=1 -q 20,20 -o output.fastq input.fastq.gz";



my $adpt_fwd = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" ;
my $adpt_rev = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"; 


# page 11 http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/scriptseq-v2-rna-seq/scriptseq-v2-rna-seq-library-prep-guide.pdf
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

# page 17 http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/scriptseq-complete/scriptseq-complete-kit-human-mouse-rat-low-input-library-prep-guide.pdf
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT



my $samtoolsBin = "samtools"; 
my $inputFolder = "./";
my $outputFolder = "./RNA-SEQ";
my $tophatProgram = "tophat";


my $version;
my $ucsc;
my $ensembl;
my $human;
my $mouse;

my %gf_hash;
my %bt2_hash;

my $gf;
my $bowind;
my $paired;
my $help;

GetOptions (
			"help"  => \$help,   # flag.
			"version=s" => \$version,    # string
            "inputFolder=s"   => \$inputFolder,      # string
            "outputFolder=s"   => \$outputFolder,      # string
            "ucsc"  => \$ucsc,   # flag. # UCSC returns gene symbols
            "ensembl"  => \$ensembl,   # flag. # Ensembl returns ensembl ids
            "human"  => \$human,   # flag.
            "mouse"  => \$mouse,   # flag.
			"paired"  => \$paired),   # flag. # default is single
or usage("Invalid commmand line options.");

if($help){ usage("See commmand line options."); }

usage("The organism must be specified.")
   unless (defined $human or defined $mouse);

usage("The namespace must be specified.")
   unless (defined $ucsc or defined $ensembl);

usage("The version must be specified.")
   unless defined $version;



sub usage {
   my $message = $_[0];
   if (defined $message && length $message) {
      $message .= "\n"
         unless $message =~ /\n$/;
   }

   my $command = $0;
   $command =~ s#^.*/##;

   print STDERR (
      $message,
      "usage: $command --inputFolder inputFolder [--outputFolder outputFolder] [--paired] --ensembl[--ucsc] --human[|--mouse] --version\n" .
      "       Organism is required\n" .
      "       Namespace is required\n" .
      "       Version of genome and annotation is required\n" .
      "       Input folder defaults to ./\n" .
      "       Output folder defaults to ./RNA-SEQ\n" .
      "       Defaults to single end. Specify --paired if needed\n"
   );

   die("\n")
}

if($ucsc){

  if($human){

	%gf_hash = (
    "18"  => "NA",
    "19" => "/nfs3/PHARM/Morgun_Lab/richrr/db/TopHat/Homo_sapiens/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf",
	);

	%bt2_hash = (
    "18"  => "NA",
    "19" => "/nfs3/PHARM/Morgun_Lab/richrr/db/TopHat/Homo_sapiens/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome",
	);
	
  } elsif($mouse){
 
	%gf_hash = (
    "9"  => "NA",
    "10" => "NA",
	);

	%bt2_hash = (
    "9"  => "NA",
    "10" => "NA",
	);
 
  }

} elsif($ensembl){

  if($human){

	%gf_hash = (
    "18"  => "NA",
    "37" => "/nfs3/PHARM/Morgun_Lab/richrr/db/TopHat/Homo_sapiens/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
	);

	%bt2_hash = (
    "18"  => "NA",
    "37" => "/nfs3/PHARM/Morgun_Lab/richrr/db/TopHat/Homo_sapiens/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome",
	);
	
  } elsif($mouse){
 
	%gf_hash = (
    "9"  => "NA",
    "37" => "/nfs3/PHARM/Morgun_Lab/darria/Genomes/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/genes.gtf",
	);

	%bt2_hash = (
    "9"  => "NA",
    "37" => "/nfs3/PHARM/Morgun_Lab/darria/Genomes/Mus_musculus/Ensembl/NCBIM37/Sequence/Bowtie2Index/genome",
	);
 
  }

}


$gf = $gf_hash{$version};
$bowind = $bt2_hash{$version};

#print $gf . "\n" . $bowind . "\n";


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	#	work flow
	#	1. remove adapters
	#	2. remove low quality sequences(better performed before remove of "N",for some file ,every sequence has an "N" in a particular possition, if we remove "N" first, all sequences will be removed. but these "N"s may be removed in this quality filter step, so after it is removed, we do not have to worry about all sequences to be discarded)
	#	3. remove sequence with any "N"
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 														remove sequence with "N"====record change of file after each step
#																	use prinseq
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

my $outFolder = "$outputFolder.AfterCutAdapt";
`mkdir $outFolder` if ( ! -d $outFolder);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


sub cutadapt{
	my ($inputFile, $OutPutFile) = @_;

    my $cutadaptTemplate = "cutadapt -a $adpt_fwd -a $adpt_rev -m 60 --max-n=1 -q 20,20 -o output.fastq input.fastq.gz";
	$cutadaptTemplate =~s/input.fastq.gz/$inputFile/g;
	$cutadaptTemplate =~s/output.fastq/$OutPutFile/g;
	print $cutadaptTemplate,"\n\n"; #die;
	`$cutadaptTemplate` ;
}


# paired end
sub cutadapt_pe{
	my ($inputFileR1, $OutPutFileR1, $inputFileR2, $OutPutFileR2) = @_;
	
	my $cutadaptTemplate = "cutadapt -a $adpt_fwd -A $adpt_rev -m 60 --max-n=1 -q 20,20 -o out.1.fastq -p out.2.fastq reads.1.fastq.gz reads.2.fastq.gz";

	$cutadaptTemplate =~s/reads.1.fastq.gz/$inputFileR1/g;
	$cutadaptTemplate =~s/out.1.fastq/$OutPutFileR1/g;

	$cutadaptTemplate =~s/reads.2.fastq.gz/$inputFileR2/g;
	$cutadaptTemplate =~s/out.2.fastq/$OutPutFileR2/g;

	print $cutadaptTemplate,"\n\n"; #die;
	`$cutadaptTemplate` ;
}

#use sequencePreprocessing::FileHandle;


sub processFastq{
	my @gzFilesR2 ;
	
	my @gzFilesR1 = `ls $inputFolder | grep "fastq.gz" | grep "lane" | grep "_R1_" `;     
    if($paired)	{
		@gzFilesR2 = `ls $inputFolder | grep "fastq.gz" | grep "lane" | grep "_R2_" `;
	}
	#print @gzFilesR1,"\n" , @gzFilesR2,"\n\n"; 
    #print "Size= ", scalar @gzFilesR1 , " ", $#gzFilesR1+1 , "\n";
        
        for(my $i=0; $i <= $#gzFilesR1; $i++){
		
			my $fileNameR2 ;
			my $InputFastqFileR2;
			my $OutPutFileR2;
			my $cutAdaptorFileR2;

	        my $fileNameR1 = $gzFilesR1[$i] ;
			chomp $fileNameR1; 
			
			if($paired)	{
				$fileNameR2 = $gzFilesR2[$i] ;        
				chomp $fileNameR2 ;
			}
			print "$fileNameR1 $fileNameR2\n"; # next;
			
			my $fileName = "$fileNameR1";
			if($paired)	{
				$fileName = "$fileNameR1-$fileNameR2";
			}
			$fileName =~s/_001.fastq.gz//g;
			
			my $InputFastqFileR1 = "$inputFolder/$fileNameR1";
			my $OutPutFileR1 = "$outFolder/$fileNameR1.cutadapt";

			if($paired)	{
				$InputFastqFileR2 = "$inputFolder/$fileNameR2";
				$OutPutFileR2 = "$outFolder/$fileNameR2.cutadapt";
				cutadapt_pe($InputFastqFileR1, $OutPutFileR1, $InputFastqFileR2, $OutPutFileR2) if ((! -e "$OutPutFileR1") && (! -e "$OutPutFileR2"));
			} else {
				# run single end cutadapt
				cutadapt($InputFastqFileR1, $OutPutFileR1) if (! -e "$OutPutFileR1");
			}

			my $cutAdaptorFileR1 = $OutPutFileR1;
			if($paired)	{
				$cutAdaptorFileR2 = $OutPutFileR2;
			}
			my $tophatOutFolder = "$outFolder.tophat/$fileName";
			`mkdir $outFolder.tophat` if (! -d "$outFolder.tophat");
			if($paired)	{
				runTopHat_pe($cutAdaptorFileR1,$cutAdaptorFileR2,$tophatOutFolder) if (! -e "$tophatOutFolder/accepted_hits.bam");
			} else {
				# run single end tophat
				runTopHat($cutAdaptorFileR1,$tophatOutFolder) if (! -e "$tophatOutFolder/accepted_hits.bam");
			}


			my $samFolder = "$outFolder.tophat.Aligned.sam";
			`mkdir $samFolder` if (! -d $samFolder);
			my $acceptedBam = "$tophatOutFolder/accepted_hits.bam";
			# sort by name, convert to SAM for htseq-count
			print "$samtoolsBin sort -n --threads 20 -o $acceptedBam.sn.bam --output-fmt BAM $acceptedBam  \n";
			`$samtoolsBin sort -n --threads 20 -o $acceptedBam.sn.bam --output-fmt bam $acceptedBam ` if (! -e "$samFolder/$fileName");
			`$samtoolsBin view --threads 20 --output-fmt sam -o $samFolder/$fileName $acceptedBam.sn.bam` if (! -e "$samFolder/$fileName");
# samtools view -?  
# 3. SAM->BAM conversion: `samtools view -b in.sam.gz'.    -o FILE  output file name
#  4. BAM->SAM conversion: `samtools view -h in.bam'.  # do you need -h?


			
			my $countFolder = "$outFolder.tophat.Aligned.sam.htseqCount";
			`mkdir $countFolder` if (! -d $countFolder);
			#my $htseqCommand = "htseq-count --order=name -a 10 -f sam --stranded=no -o $countFolder/$fileName $samFolder/$fileName $gf"; # outputs SAM file
			my $htseqCommand = "htseq-count --order=name -a 10 -f sam --stranded=no $samFolder/$fileName $gf > $countFolder/$fileName.txt"; 
			
			print "\n\n$htseqCommand\n\n\n"; #die;
			`$htseqCommand\n` if (! -e "$countFolder/$fileName");


		} 
}

sub runTopHat_pe {
	my ($cutAdaptorFileR1, $cutAdaptorFileR2, $tophatOutFolder) = @_;
	my $cmd = "$tophatProgram --no-coverage-search -G $gf -p 20 -o $tophatOutFolder $bowind $cutAdaptorFileR1 $cutAdaptorFileR2 ";  
	print ($cmd,"\n\n");
	`$cmd`;	
}


sub runTopHat {
	my ($cutAdaptorFile, $tophatOutFolder) = @_;
	my $cmd = "$tophatProgram --no-coverage-search -G $gf -p 8 -o $tophatOutFolder $bowind $cutAdaptorFile";
	print ($cmd,"\n\n");
	`$cmd`;	
}
processFastq;

