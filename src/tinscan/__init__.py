#!/usr/bin/env python
from collections import namedtuple

from .LASTZ_wrapper import *


def readLASTZ(infile, minID=90):
    """Read in LASTZ result file from LASTZ_genome_align.sh
    Populate nested dictionary of hits keyed by Target and then Query scaffold names."""
    with open(infile) as f:
        content = f.readlines()
    content = [x.strip().split() for x in content]
    hitsDict = dict()
    # Set named tuple format
    hitTup = namedtuple(
        "Elem",
        [
            "t_start",
            "t_end",
            "t_strand",
            "q_start",
            "q_end",
            "q_strand",
            "idPct",
            "UID",
        ],
    )
    counter = 0
    # Read in split rows
    for row in content:
        # Ignore lines begining with '#'
        if row[0][0] == "#":
            continue
        # Count if hit identity exceeds threshold
        elif float(row[9]) >= minID:
            counter += 1
            UID = counter
            t_name = str(row[0])
            t_strand = str(row[1])
            t_start = int(row[2]) - 1  # Convert from idx '1' to idx '0'
            t_end = int(row[3]) - 1  # Convert from idx '1' to idx '0'
            q_name = str(row[4])
            q_strand = str(row[5])
            idPct = float(row[9])
            # Check that start position < end position
            if int(row[6]) < int(row[7]):
                q_start = int(row[6]) - 1  # Convert from idx '1' to idx '0'
                q_end = int(row[7]) - 1  # Convert from idx '1' to idx '0'
            # Correct for inverted coordinates
            else:
                print("Inverting query sequence coordinates for record number: ", UID)
                q_end = int(row[6]) - 1  # Convert from idx '1' to idx '0'
                q_start = int(row[7]) - 1  # Convert from idx '1' to idx '0'
            # Create target scaffold dict if not seen
            if t_name not in hitsDict.keys():
                hitsDict[t_name] = dict()
            # Create query scaffold list if not seen
            if q_name not in hitsDict[t_name].keys():
                hitsDict[t_name][q_name] = list()
            # Write record to target:query list as named tuple
            hitsDict[t_name][q_name].append(
                hitTup(t_start, t_end, t_strand, q_start, q_end, q_strand, idPct, UID)
            )
    return hitsDict


def getHitPairs(hits, args):
    """Given a nested dictionary keyed by Target scaffold name, then Query scaffold name,
    where Query scaffold sub-dict contains a list of hits stored as named tuples
    i.e. (t_start,t_end,t_strand,q_start,q_end,q_strand,idPct,UID)
    return pairs of hits which satisfy the filter criteria as a list of objects with structure:
    (
    ("Target_name","Query_name"),
    hit=hitTup(t_start,t_end,t_strand,q_start,q_end,q_strand,idPct,UID),
    mate=hitTup(t_start,t_end,t_strand,q_start,q_end,q_strand,idPct,UID)
    )
    - that is a tuple containg three tuples."""
    valid_elem = []
    for t_name in hits.keys():
        for q_name in hits[t_name].keys():
            for hit in hits[t_name][q_name]:
                for mate in hits[t_name][q_name]:
                    if (
                        mate.UID != hit.UID
                        and mate.q_strand == hit.q_strand
                        and abs(mate.idPct - hit.idPct) <= args.maxIdentDiff
                        and mate.t_start - hit.t_end >= args.minInsert
                        and mate.t_start - hit.t_end <= args.maxInsert
                    ):
                        if (
                            hit.q_strand == "+"
                            and mate.q_start - hit.q_end <= args.qGap
                            and mate.q_start - hit.q_end >= 0 - args.maxTSD
                        ):
                            valid_elem.append(((t_name, q_name), hit, mate))
                        elif (
                            hit.q_strand == "-"
                            and hit.q_start - mate.q_end <= args.qGap
                            and hit.q_start - mate.q_end >= 0 - args.maxTSD
                        ):
                            valid_elem.append(((t_name, q_name), hit, mate))
    return valid_elem


def formatGFFline(pair, featureID):
    seqid = str(pair[0][0])
    source = "InsertScanner"
    feature_type = "Candidate_Insertion"
    start = pair[1].t_end + 1
    end = pair[2].t_start - 1
    score = "."
    strand = pair[1].t_strand
    phase = "."
    leftflank = (
        pair[0][1]
        + "_"
        + '"'
        + pair[1].q_strand
        + '"'
        + "_"
        + str(pair[1].q_start)
        + "_"
        + str(pair[1].q_end)
    )
    rightflank = (
        pair[0][1]
        + "_"
        + '"'
        + pair[2].q_strand
        + '"'
        + "_"
        + str(pair[2].q_start)
        + "_"
        + str(pair[2].q_end)
    )
    attributes = (
        "ID="
        + str(featureID)
        + ";len="
        + str(end - start)
        + ";leftflank="
        + leftflank
        + ";rightflank="
        + rightflank
        + ";leftID="
        + str(pair[1].idPct)
        + ";rightID="
        + str(pair[2].idPct)
        + "\n"
    )
    return "\t".join(
        [
            seqid,
            source,
            feature_type,
            str(start),
            str(end),
            score,
            strand,
            phase,
            attributes,
        ]
    )


def getTSD(pair, parentID):
    TSDlen = None
    if pair[1].q_strand == "+" and pair[1].q_end > pair[2].q_start:
        TSDlen = pair[1].q_end - pair[2].q_start
    elif pair[1].q_strand == "-" and pair[1].q_start < pair[2].q_end:
        TSDlen = pair[2].q_end - pair[1].q_start
    if TSDlen:
        TSDr_start = pair[2].t_start
        TSDr_end = pair[2].t_start + TSDlen
        TSDl_start = pair[1].t_end - TSDlen
        TSDl_end = pair[1].t_end
        strand = pair[1].t_strand
        seqid = str(pair[0][0])
        source = "InsertScanner"
        feature_type = "Candidate_TSD"
        attributesR = (
            "ID="
            + str(parentID)
            + "_TSD_R"
            + ";Parent="
            + str(parentID)
            + ";len="
            + str(TSDlen)
            + "\n"
        )
        attributesL = (
            "ID="
            + str(parentID)
            + "_TSD_L"
            + ";Parent="
            + str(parentID)
            + ";len="
            + str(TSDlen)
            + "\n"
        )
        # Yield left
        yield "\t".join(
            [
                seqid,
                source,
                feature_type,
                str(TSDl_start),
                str(TSDl_end),
                ".",
                strand,
                ".",
                attributesL,
            ]
        )
        # Yield right
        yield "\t".join(
            [
                seqid,
                source,
                feature_type,
                str(TSDr_start),
                str(TSDr_end),
                ".",
                strand,
                ".",
                attributesR,
            ]
        )
    else:
        return None


def getFlank(pair, parentID):
    """pair structure is:
    pair[0] = ("Target_name","Query_name")
    pair[1] = hit = left = hitTup(t_start,t_end,t_strand,q_start,q_end,q_strand,idPct,UID)
    pair[2] = mate = right = hitTup(t_start,t_end,t_strand,q_start,q_end,q_strand,idPct,UID)
    """
    seqid = str(pair[0][0])
    source = "InsertScanner"
    feature_type = "InsertFlank"
    score = "."
    phase = "."
    # Compose Leftflank
    strandL = pair[1].t_strand
    leftflank = (
        pair[0][1]
        + "_"
        + '"'
        + pair[1].q_strand
        + '"'
        + "_"
        + str(pair[1].q_start)
        + "_"
        + str(pair[1].q_end)
    )
    attributesL = (
        "ID="
        + str(parentID)
        + "_Flank_L"
        + ";Parent="
        + str(parentID)
        + ";len="
        + str(pair[1].t_end - pair[1].t_start)
        + ";leftflank="
        + leftflank
        + ";leftID="
        + str(pair[1].idPct)
        + "\n"
    )
    # Compose RightFlank
    strandR = pair[2].t_strand
    rightflank = (
        pair[0][1]
        + "_"
        + '"'
        + pair[2].q_strand
        + '"'
        + "_"
        + str(pair[2].q_start)
        + "_"
        + str(pair[2].q_end)
    )
    attributesR = (
        "ID="
        + str(parentID)
        + "_Flank_R"
        + ";Parent="
        + str(parentID)
        + ";len="
        + str(pair[2].t_end - pair[2].t_start)
        + ";rightflank="
        + rightflank
        + ";rightID="
        + str(pair[2].idPct)
        + "\n"
    )
    # Yield LeftFlank
    yield "\t".join(
        [
            seqid,
            source,
            feature_type,
            str(pair[1].t_start),
            str(pair[1].t_end),
            score,
            strandL,
            phase,
            attributesL,
        ]
    )
    # Yield RightFlank
    yield "\t".join(
        [
            seqid,
            source,
            feature_type,
            str(pair[2].t_start),
            str(pair[2].t_end),
            score,
            strandR,
            phase,
            attributesR,
        ]
    )


def writeGFFlines(validPairs, reportFlanks=True):
    yield "#gff-version 3\n#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes\n"
    counter = 0
    fillLen = len(str(abs(len(validPairs)))) + 1
    for pair in validPairs:
        counter += 1
        featureID = "IS_" + str(counter).zfill(fillLen)
        yield formatGFFline(pair, featureID)
        TSDlines = getTSD(pair, featureID)
        if TSDlines:
            for x in TSDlines:
                yield x
        if reportFlanks:
            for y in getFlank(pair, featureID):
                yield y
