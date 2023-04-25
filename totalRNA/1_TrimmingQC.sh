##########################################################
### Quality control and trimming of reads (paired end) ###
##########################################################

declare -r PROGRAMDIR="/share/tools/"
declare SOURCEDIR="/ngs-data/data/hecatos/Hepatic/Azathioprine/TotalRNA/concatenated/"   #location of raw data
declare SEQMODE="paired"   #specify "paired", "single"
declare SUFFIX_in='.fastq.gz' #specify extension of input files (".fastq" or ".fastq.gz")

#For paired end mode: specify the SUFFIXES used for read1, read2 (optional) 
declare PAIR1="_R1"
declare PAIR2="_R2"

#############################################
### Defining parameters for script to run ### 
#############################################

declare OUTPUTDIR="/share/analysis/Carlo/trimming/"	# Where all files will be output
declare TRIMM_DIR="${OUTPUTDIR}1_Trimmed/"
declare QC_fastp="${OUTPUTDIR}2_fastpQC_output/"
declare multiQC_fastp="${OUTPUTDIR}2_MultiQC_fastp/"

declare SUFFIX_out="_trimmed${SUFFIX_in}"

mkdir -p ${OUTPUTDIR}
mkdir -p ${TRIMM_DIR}
mkdir -p ${QC_fastp}
mkdir -p ${multiQC_fastp}

#These Echos are for me ensuring the parameters are set properly
echo "Reading files from $SOURCEDIR"
echo "Outputting files to $OUTPUTDIR"


##################################
### 1 - Fastp : TRIM RAW READS ###
##################################

# trimming 
# fastp v0.19.8

declare FILES1=${SOURCEDIR}*${PAIR1}${SUFFIX_in}

for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME}	
	READ2=${FILENAME:0:-${#PAIR1}-${#SUFFIX_in}}${PAIR2}${SUFFIX_in}	
	echo ${FILENAME}
	echo -e "[TRIMMING] fastp: [${READ1:${#SOURCEDIR}:-${#PAIR1}-${#SUFFIX_in}}]" 
	fastp \
	--in1 ${READ1} --in2 ${READ2} \
	--out1 ${TRIMM_DIR}${READ1:${#SOURCEDIR}:-${#SUFFIX_in}}${SUFFIX_out} \
	--out2 ${TRIMM_DIR}${READ2:${#SOURCEDIR}:-${#SUFFIX_in}}${SUFFIX_out} \
	--json ${QC_fastp}${READ1:${#SOURCEDIR}:-${#SUFFIX_in}-${#PAIR1}}"PE_fastp.json" \
	--html ${QC_fastp}${READ1:${#SOURCEDIR}:-${#SUFFIX_in}-${#PAIR1}}"PE_fastp.html" \
	--cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 30 --length_required 36 
done
"/share/analysis/Carlo/trimming/QC.fastp/"
echo "Trimming done!"

################################
#### If additional trimming is needed (see multiQCreport): 
# add to R1: --trim_front1 {amount_bases} and/or --trim_tail1 {amount_bases}
# add to R2: --trim_front2 {amount_bases} and/or --trim_tail2 {amount_bases}


#################################################
### Quality control raw reads: MultiQC report ###
#################################################

# Running multiQC on fastp output

multiqc ${QC_fastp} --filename MultiQC_Report.html --outdir ${multiQC_fastp}

echo "MultiQC report for trimming created!"
