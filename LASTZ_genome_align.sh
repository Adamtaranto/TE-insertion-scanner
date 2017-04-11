#!/bin/sh

# Reset OPT index to 1
OPTIND=1         

# Initialize our own variables:
minIdt="60"
minLen="100"
outdir="outdir"
prefix="genome_alignment"
LZ="lastz"

while getopts ":d:p:t:q:i:l:z:" opt; do #Options followed by ":" expect an arguement.
	case "$opt" in
	z)	LZ=$OPTARG
		if [ ! -f "$LZ" ]; then
			echo "LASTZ does not exist at path: $LZ" >&2
			exit 1
		fi
		echo "Path to LASTZ app set as: $OPTARG" >&2
		;;
	i)	minIdt=$OPTARG
		echo "Setting min identity threshold at $OPTARG %" >&2
		;;
	l)	minLen=OPTARG
		echo "Setting min hit length threshold to $OPTARG" >&2
		;;
	d)	outdir=$OPTARG
		;;
	p)	prefix=$OPTARG
		;;
	t)	tdata=$OPTARG
		if [ ! -d "$tdata" ]; then
			echo "Target sequence directory $tdata does not exist." >&2
			exit 1
		fi
		;;
	q)	qdata=$OPTARG
		if [ ! -d "$qdata" ]; then
			echo "Query sequence directory $qdata does not exist." >&2
			exit 1
		fi
		;;
	\?)	echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "output_file='$output_file', Leftovers: $@"

# Housekeeping
if [ ! -d "$outdir" ]; then
	mkdir $outdir # Make output directory if does not exist
fi

# Make empty output files
out_gff=$(echo $outdir"/"$prefix".gff3")
out_tab=$(echo $outdir"/"$prefix"_concat.tab")

# Scrub old results
if [ -f "$out_gff" ]; then
	echo "Removing old output file: $out_gff"
	rm $out_gff
fi

# Scrub old results
if [ -f "$out_tab" ]; then
	echo "Removing old output file: $out_tab"
	rm $out_tab
fi

# Initialise output files
echo -e "##gff-version 3\n#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes\n" > $out_gff
echo -e "#name1\tstrand1\tstart1\tend1\tname2\tstrand2\tstart2+\tend2+\tscore\tidentity\n" > $out_tab

# Run alignments
for target in $(find $tdata -regex ".*\.\(fa\|fasta\|fsa\)" | sort)
do
	for query in $(find $qdata -regex ".*\.\(fa\|fasta\|fsa\)" | sort)
	do
	t_file=$(basename $target)
	t_name="${t_file%.*}"
	q_file=$(basename $query)
	q_name="${q_file%.*}"
	outfile=$(echo $outdir"/"$t_name".vs."$q_name".tab")
	outfile_filtered=$(echo $outfile"_filtered.bed")
	outfile_filtered_sorted=$(echo $outfile"_filtered_sorted.bed")
	feature=$(echo $q_name"_ident")
	$LZ $target $query \
	--gfextend \
	--chain \
	--gapped \
	--step=1 \
	--strand=both \
	--hspthresh=3000 \
	--entropy \
	--output=$outfile \
	--format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,score,identity \
	--verbosity=1 \
	--markend
	#Scrub % symbols
	sed -i '' -e 's/%//g' $outfile
	## Filter Inter_Chrome targets to min len $minLen [100], min identity $minIdt [90]
	## New fields = name1,strand1,start1,end1,name2,strand2,start2+,end2+,score,identity
	awk '!/^#/ { print; }' $outfile | awk -v minLen="$minLen" '0+$5 >= minLen {print ;}' | awk -v OFS='\t' -v minIdt="$minIdt" '0+$13 >= minIdt {print $1,$2,$3,$4,$6,$7,$8,$9,$11,$13;}' | sed 's/ //g' > $outfile_filtered
	## Sort filtered bed file by chrom, start, stop
	sort -k 1,1 -k 3n,4n $outfile_filtered >> $out_tab
	rm $outfile
	rm $outfile_filtered
	done
done

## Create GFF3 file for merged filtered hits
awk -v OFS='\t' '!/^#/ BEGIN{i=0}{i++;}{j= sprintf("%09d", i)}{print $1,"LASTZ","'"$feature"'",$3,$4,$9,$2,".","ID=LZ_Targets_"j";Idt="$10";Target="$5"_"$6"_"$7"_"$8 ;}' $out_tab >> $out_gff

# End of file