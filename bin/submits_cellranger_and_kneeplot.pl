#!/usr/bin/perl

############################
#### Step 1 -- Define parameters, default modules, infiles, and other dependencies
############################

use LoadParameters::Parameters;
use ReformatPerlEntities::ObtainOutfileWOpath;
use PathsDefinition::PathsToInputs;

$ThisProgramName = $0;
$ThisProgramName =~ s/\S+\///;

$CommentsForHelp = "
#####################################################################################
################### START INSTRUCTIONS TO RUN THIS PROGRAM ##########################
###
### Entering *fastq files runs 'cellranger count', then obtains a knee-plot
### analysis to get a list of cell barcodes alternative to the cellranger approach,
### used knee-plot based list of cell barcodes and runs 'cellranger reanalyze' with them
###
### -------------------------------------------INPUTS----------------------------------------------
###
### [1]
### -indir_fastq indicates the directory containing *fastq files
###
### ----------------------------------------MAIN OUTPUTS-------------------------------------------
###
### [1]
### All files from a 'cellranger count' run, including:
### .../SC_RNA_COUNTER_CS
### .../outs/analysis
### .../outs/filtered_gene_bc_matrices
### .../outs/raw_gene_bc_matrices
### .../outs/kneeplot
### .../outs/selected_gene_bc_matrices
### .../outs/reanalysis
### and log files
###
### ------------------------------------------COMMANDS---------------------------------------------
###
### $ThisProgramName [options]
###   -path_outfiles      (path/name to the directory where outfiles will be saved)
###   -indir_fastq        (path to the directory where fastq files live, like:
###                         '/path_to/112233_A00123_0123_ABCD1EFGXX_Javier' in files
###                         '/path_to/112233_A00123_0123_ABCD1EFGXX_Javier/Sample_1_S1_L001_R1_001_fastq.gz'
###                         '/path_to/112233_A00123_0123_ABCD1EFGXX_Javier/Sample_1_S1_L001_R2_001_fastq.gz'
###
###   -prefix_fastq        (prefix the fastq files before the 'SX_L00Y' sample/lane numbers, like:
###                         'Sample_1' in file names above
###                         this is useful when multiple sample fastq files are contained into the same -indir_fastq
###
###   -sample_id           (prefix for the outfiles)
###
###   -numb_cells          (number of expected cells for --expect-cells for cellranger (e.g. '3000')
###                         Note independently of this parameter, a knee-pot will be used to resubmit to cellranger)
###
###   -step_numbers_to_run (<comma> and <dash> delimited number of steps to run, e.g. '1,2' or '3' or 'ALL')
###                          Steps
###                          1: Runs cellranger
###                          2: Gets number of cells from knee-plot
###                          3: Runs re-analysis by cellranger
###
###   -cluster_name        (name of either of the following clusters where process are running)
###                         'mordor'
###                         'h4h'
###                         'samwise'
###                         'local'
###
###   -genome_reference    (number of either of the following genome references to map reads)
###                         '1' for 'GRCh38-1.2.0' (human)
###                         '2' for 'GRCh38-1.2.0_premrna' (human Nuc-Seq)
###                         '3' for 'mm10-1.2.0' (mouse)
###                         '4' for 'mm10-1.2.0_premrna' (mouse Nuc-Seq)
###
##################### END INSTRUCTIONS TO RUN THIS PROGRAM ##########################
#####################################################################################";
$CommentsForHelp =~ s/(\n####)|(\n###)|(\n##)|(\n#)/\n/g;

&Readme;
&Parameters;

$maxjobs = 40; ## maximum number of steps to submit at once, each step will take a core
$mempercore = 10; ## requested memory per core

if ($hashParameters{numb_cells} =~ /\d/) {
$cells_for_knee = $hashParameters{numb_cells} * 10;
}else{
die "\n\nERROR!!! unexpected format of -numb_cells $hashParameters{numb_cells}\n\n";
}


############################
#### Step 2 -- Define and check dependencies
############################

#### Note: dependencies must be visible from the node submitting the jobs

### Software
$dropseq_bth        = "~/PROGRAMS/DROPSEQ/Drop-seq_tools-1.13/BAMTagHistogram";
$cellranger_exe     = "~/PROGRAMS/CELLRANGER/cellranger-2.1.1/cellranger";
$r_exe              = "R"; # R will be loaded by module load R

### References
$References{1}      = "~/DATABASES/10X/refdata-cellranger-GRCh38-1.2.0";
$References{2}      = "~/DATABASES/10X/refdata-cellranger-GRCh38-1.2.0_premrna";
$References{3}      = "~/DATABASES/10X/refdata-cellranger-mm10-1.2.0";
$References{4}      = "~/DATABASES/10X/refdata-cellranger-mm10-1.2.0_premrna";
$transcriptome = $hashParameters{$hashParameters{genome_reference}};

### Check that dependencies are available
$dropseq_bth    =~ s/~\//\/$Users_home\/$DefaultUserName\//;
$cellranger_exe =~ s/~\//\/$Users_home\/$DefaultUserName\//;
$transcriptome  =~ s/~\//\/$Users_home\/$DefaultUserName\//;

unless (-f $dropseq_bth) {
die "\n\nERROR!!! couldn't find '$dropseq_bth' (BAMTagHistogram)";
}
unless (-f $cellranger_exe) {
die "\n\nERROR!!! couldn't find '$cellranger_exe' (cellranger)";
}
unless (-d $transcriptome) {
die "\n\nERROR!!! couldn't find '$transcriptome' (transcriptome)";
}

### Get job mode
if ($hashParameters{cluster_name} =~ /^mordor$/) {
$jobmode      = "sge";
$queue_params = "qsub -q highmem.q";
}elsif ($hashParameters{cluster_name} =~ /^h4h$/) {
$jobmode = "torque";
$queue_params = "qsub -q himem -l vmem=60G,walltime=72:00:00";
}elsif ($hashParameters{cluster_name} =~ /^samwise$/) {
$jobmode = "local";
$queue_params = "";
}else{
die "\n\nERROR!!! unrecognized -cluster_name '$hashParameters{cluster_name}'\n\n";
}

############################
#### Step 3 -- Get list of steps to run
############################
if ($hashParameters{step_numbers_to_run} =~ /^ALL$/i) {
	foreach $s (1..3) {
	$hash_step_numbers_to_run{$s} = 1;
	}
}

############################
#### Step 4 -- Send to cellranger
############################

if ($hash_step_numbers_to_run{1}) {
open  INSFOR_CELLRANGER, ">$hashParameters{path_outfiles}/$hashParameters{sample_id}.run_cellranger.sh" or die "Can't open '$hashParameters{path_outfiles}/$hashParameters{sample_id}.run_cellranger.sh'\n";
print INSFOR_CELLRANGER "
#!/bin/bash
#
cd $hashParameters{path_outfiles}
#
$cellranger_exe count \
--id=$hashParameters{prefix_fastq} \
--sample=$hashParameters{sample_id} \
--transcriptome=$transcriptome \
--fastqs=$hashParameters{indir_fastq} \
--expect-cells=$hashParameters{numb_cells} \
--jobmode=$jobmode \
--mempercore=$mempercore \
--maxjobs=$maxjobs\n";
close INSFOR_CELLRANGER;

### Runs the script
system "chmod +x $hashParameters{path_outfiles}/$hashParameters{sample_id}.run_cellranger.ins";
system "$queue_params $hashParameters{path_outfiles}/$hashParameters{sample_id}.run_cellranger.ins";
}


############################
#### Step 5 -- Get number of reads per barcode and knee plot from cellranger's possorted_genome_bam.bam
############################

if ($hash_step_numbers_to_run{2}) {

### Instructions to get number of reads per barcode and run knee plot script
$expected_possorted_genome_bam = "$hashParameters{path_outfiles}/possorted_genome_bam.bam";
$expected_possorted_genome_bam =~ s/~\//\/$Users_home\/$DefaultUserName\//;
$possorted_genome_reads_per_barcode = "$hashParameters{path_outfiles}/possorted_genome_reads_per_barcode.tsv";
$possorted_genome_reads_per_barcode =~ s/~\//\/$Users_home\/$DefaultUserName\//;

	unless (-f $expected_possorted_genome_bam) {
	die "\nERROR!!! couldn't find '$expected_possorted_genome_bam'\n\n";
	}
	
open  INSFOR_KNEEPLOT, ">$hashParameters{path_outfiles}/$hashParameters{sample_id}.run_bamtaghistogram.sh" or die "Can't open '$hashParameters{path_outfiles}/$hashParameters{sample_id}.run_bamtaghistogram.sh'\n";
print INSFOR_KNEEPLOT "module load java/8
cd $hashParameters{path_outfiles}
$dropseq_bth TAG=CB I=$expected_possorted_genome_bam O=$possorted_genome_reads_per_barcode
#
module load R
$r_exe --no-save < $hashParameters{path_outfiles}/$hashParameters{sample_id}.run_dropbead.R
mv ~/$hashParameters{sample_id}_kneeplot.pdf $hashParameters{path_outfiles}\n";
close INSFOR_KNEEPLOT;

### Instructions to  make knee-plot
open  INSFOR_DROPDEAD, ">$hashParameters{path_outfiles}/$hashParameters{sample_id}.run_dropbead.R" or die "Can't open '$hashParameters{path_outfiles}/$hashParameters{sample_id}.run_dropbead.R'\n";
print INSFOR_DROPDEAD "library(dropbead)
mat<-as.data.frame(read.table(\"$possorted_genome_reads_per_barcode\"),header=F)
pdf(\"~/$hashParameters{sample_id}_kneeplot.pdf\")
plotCumulativeFractionOfReads(mat, cutoff = $cells_for_knee, draw.knee.point = TRUE)
InflectionPoint<-estimateCellNumber(mat[, 1], max.cells = $cells_for_knee)
write(file=\"$hashParameters{path_outfiles}/$hashParameters{sample_id}_kneeplot_inflection.txt\",x=InflectionPoint)
dev.off()
q()\n";
close INSFOR_DROPDEAD;
	
### Runs the script
system "chmod +x $hashParameters{path_outfiles}/$hashParameters{sample_id}.run_bamtaghistogram.ins";
system "chmod +x $hashParameters{path_outfiles}/$hashParameters{sample_id}.run_dropbead.ins";
system "$queue_params $hashParameters{path_outfiles}/$hashParameters{sample_id}.run_bamtaghistogram.ins";

}


############################
#### Step 6 -- Get the list of barcode ID's based on the knee plot inflection point
############################

if ($hash_step_numbers_to_run{3}) {
$expected_knee_inflection = "$hashParameters{path_outfiles}/possorted_genome_bam.bam$hashParameters{path_outfiles}/$hashParameters{sample_id}_kneeplot_inflection.txt";
$expected_knee_inflection =~ s/~\//\/$Users_home\/$DefaultUserName\//;
	
	unless (-f $expected_knee_inflection) {
	die "\nERROR!!! couldn't find '$expected_knee_inflection'\n\n";
	}

open  INSFOR_RERUN_CELLRANGER, ">$hashParameters{path_outfiles}/$hashParameters{sample_id}.run_bamtaghistogram.sh" or die "Can't open '$hashParameters{path_outfiles}/$hashParameters{sample_id}.run_bamtaghistogram.sh'\n";
print INSFOR_KNEEPLOT "module load java/8
cd $hashParameters{path_outfiles}
$dropseq_bth TAG=CB I=$expected_possorted_genome_bam O=$possorted_genome_reads_per_barcode
#
module load R
$r_exe --no-save < $hashParameters{path_outfiles}/$hashParameters{sample_id}.run_dropbead.R
mv ~/$hashParameters{sample_id}_kneeplot.pdf $hashParameters{path_outfiles}\n";
close INSFOR_KNEEPLOT;



CONTINUE HERE





&PrintParameters;

print "\n\n  Done!!!\n  Check '$hashParameters{path_outfiles}/$outfileWOpath.*' for outfiles\n\n";

exit;

########################################################
################ END OF PROGRAM ########################
########################################################


########################################################
################ START SUBROUTINES #####################
########################################################


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub Parameters {

########## Print "Usage" for user

print "$CommentsForHelp\n\n";

##########################
######## Options and Infiles

use Cwd 'abs_path';
$ScriptName = abs_path($0);
$Parameters .= "$ScriptName\n";

chomp @ARGV;
@arrayInputtedOneLineCommands = @ARGV;

%hashParametersTolookFor = (
'path_outfiles' => 1,
'infile' => 1,
);

#######################
#### Starts -- Evaluate parameters

&LoadParameters::Parameters::MainSubParameters(\%hashParametersTolookFor,\@arrayInputtedOneLineCommands);
$Parameters .= "$MoreParameters";

## Defining prefix string for OUTFILE
ReformatPerlEntities::ObtainOutfileWOpath::ObtainOutfileWOpath($hashParameters{infile});

#### Ends -- Evaluate parameters
#######################

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub PrintParameters {

### Printing out parameters. Need to be concatenated at sub(Parameters)
open PARAMETERS, ">$hashParameters{path_outfiles}/$outfileWOpath.ZZZ.Parameters" or die "Can't open '$hashParameters{path_outfiles}/$outfileWOpath.ZZZ.Parameters'\n";
print PARAMETERS "$Parameters";
close PARAMETERS;
}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sub Readme {

my ($date) = `date`;
chomp $date;

$Parameters .= "
#################################################################
# Javier Diaz -- $date
# javier.diazmejia\@gmail.com
#################################################################\n
$Extras
#################################################################
######################### PARAMETERS ############################
#################################################################\n\n";

}
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
