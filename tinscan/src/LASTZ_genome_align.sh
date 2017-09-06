#!/bin/sh

# Reset OPT index to 1
OPTIND=1         

# Initialize our own variables:
minIdt="60"
minLen="100"
outdir="outdir"
prefix="genome_alignment"
LZ="/usr/local/bin/lastz"
recycleTAB=0

while getopts "d:p:t:q:i:l:z:r" opt; do #Options followed by ":" expect an arguement.
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
	l)	minLen=$OPTARG
		echo "Setting min hit length threshold to $OPTARG" >&2
		;;
	d)	outdir=$OPTARG
		echo "Setting output directory to $OPTARG" >&2
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
	r)	recycleTAB=1 
		echo "Use existing alignment data if found." >&2
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

echo "Leftover arguements:" "$@" >&2

# Housekeeping
if [ ! -d "$outdir" ]; then
	echo "Creating output directory $outdir" >&2
	mkdir "$outdir" # Make output directory if does not exist
fi

# Make output file paths
out_gff=$(echo "$outdir""/""$prefix""_raw_alignments.gff3")
out_tab=$(echo "$outdir""/""$prefix""_concat.tab")

# Scrub old results
if [ -f "$out_tab" ] && [ $recycleTAB != 1 ]; then
	echo "Removing old output file: $out_tab" >&2
	rm "$out_tab"
	# Initialise output files
	echo $'#name1\tstrand1\tstart1\tend1\tname2\tstrand2\tstart2+\tend2+\tscore\tidentity' > "$out_tab"
elif [ ! -f "$out_tab" ]; then
	# Initialise output files
	echo $'#name1\tstrand1\tstart1\tend1\tname2\tstrand2\tstart2+\tend2+\tscore\tidentity' > "$out_tab"
	recycleTAB=0
else
	echo "Alignment file alreay exists, no further action: $out_tab" >&2
fi

# Scrub old raw hits GFF if new alignments will be generated 
if [ -f "$out_gff" ] && [ $recycleTAB != 1 ]; then #recycle == 0 == make new aligns + write new gff
	echo "Removing old output file: $out_gff" >&2
	rm "$out_gff"
fi

if [ $recycleTAB != 1 ]; then
	# Initialise raw alignment GFF
	echo $'##gff-version 3\n#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes' > "$out_gff"
	# Run alignments
	for target in $(find "$tdata" -type f \( -name "*.fa" -or -name "*.fasta" -or -name "*.fsa" \) | sort)
	do
		for query in $(find "$qdata" -type f \( -name "*.fa" -or -name "*.fasta" -or -name "*.fsa" \) | sort)
		do
		t_file=$(basename "$target")
		t_name="${t_file%.*}"
		q_file=$(basename "$query")
		q_name="${q_file%.*}"
		outfile=$(echo "$outdir""/""$t_name"".vs.""$q_name"".tab")
		$LZ "$target" "$query" \
		--gfextend \
		--chain \
		--gapped \
		--step=1 \
		--strand=both \
		--hspthresh=3000 \
		--entropy \
		--output="$outfile" \
		--format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,score,identity \
		--verbosity=1 \
		--markend
		#Scrub % symbols
		sed -i '' -e 's/%//g' "$outfile"
		## Filter Inter_Chrome targets to min len $minLen [100], min identity $minIdt [90]
		## New fields = name1,strand1,start1,end1,name2,strand2,start2+,end2+,score,identity
		## Sort filtered bed file by chrom, start, stop
		echo "Alignment finished, writing hits: $out_tab" >&2
		awk '!/^#/ { print; }' "$outfile" | awk -v minLen="$minLen" '0+$5 >= minLen {print ;}' | awk -v OFS='\t' -v minIdt="$minIdt" '0+$13 >= minIdt {print $1,$2,$3,$4,$6,$7,$8,$9,$11,$13;}' | sed 's/ //g' | sort -k 1,1 -k 3n,4n >> "$out_tab"
		## Create GFF3 file for merged filtered hits
		echo "Writing filtered hits to gff3: $out_gff" >&2
		awk '!/^#/ {print ;}' "$out_tab" | awk -v OFS='\t' -v q_name="$q_name" 'BEGIN{i=0}{i++;}{j=sprintf("%09d",i)}{print $1,"LASTZ","Raw_Alignment",$3,$4,$9,$2,".","ID=LZ_hit_"q_name"_"j";Idt="$10";Target="$5"_"$6"_"$7"_"$8 ;}' >> "$out_gff"
		rm "$outfile"
		done
	done
fi
# End of file