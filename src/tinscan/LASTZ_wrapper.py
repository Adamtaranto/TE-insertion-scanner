from shlex import quote
import glob
import os
import shutil
import subprocess
import sys
import tempfile


class Error(Exception):
    pass


def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s


def import_pairs(file=None, Adir=None, Bdir=None):
    pairs = list()
    with open(file) as f:
        for line in f:
            li = line.strip()
            if not li.startswith("#"):
                A, B = li.split()[:2]
                pairs.append((os.path.join(Adir, A), os.path.join(Bdir, B)))
    return pairs


def get_all_pairs(Adir=None, Bdir=None):
    pairs = list()
    for A in glob.glob(os.path.join(Adir, "*")):
        for B in glob.glob(os.path.join(Bdir, "*")):
            pairs.append((A, B))
    return pairs


def LASTZ_cmds(
    lzpath="lastz",
    pairs=None,
    minIdt=60,
    minLen=100,
    hspthresh=3000,
    outfile=None,
    verbose=False,
):
    if verbose:
        verb = 1
    else:
        verb = 0
    cmds = list()
    # Write header
    cmds.append(
        " ".join(
            [
                "echo $'#name1\\tstrand1\\tstart1\\tend1\\tname2\\tstrand2\\tstart2+\\tend2+\\tscore\\tidentity' >",
                quote(outfile),
            ]
        )
    )
    for A, B in pairs:
        t_file = A
        t_name = os.path.splitext(os.path.basename(A))[0]
        q_file = B
        q_name = os.path.splitext(os.path.basename(B))[0]
        temp_outfile = "_".join(["temp", q_name, "onto", t_name, ".tab"])
        # Compose LASTZ command
        cmds.append(
            " ".join(
                [
                    quote(lzpath),
                    quote(t_file),
                    quote(q_file),
                    "--entropy --format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,score,identity --markend --gfextend --chain --gapped --step=1 --strand=both --hspthresh="
                    + str(hspthresh),
                    "--output=" + temp_outfile,
                    "--verbosity=" + str(verb),
                ]
            )
        )
        # Scrub % symbols
        cmds.append(" ".join(["sed -i '' -e 's/%//g'", temp_outfile]))
        ## Filter Inter_Chrome targets to min len $minLen [100], min identity $minIdt [90]
        ## New Header = name1,strand1,start1,end1,name2,strand2,start2+,end2+,score,identity
        ## Sort filtered file by chrom, start, stop
        cmds.append(
            " ".join(
                [
                    "awk '!/^#/ { print; }'",
                    temp_outfile,
                    "| awk -v minLen=" + str(minLen),
                    "'0+$5 >= minLen {print ;}' | awk -v OFS=" + "'\\t'",
                    "-v minIdt=" + str(minIdt),
                    "'0+$13 >= minIdt {print $1,$2,$3,$4,$6,$7,$8,$9,$11,$13;}' | sed 's/ //g' | sort -k 1,1 -k 3n,4n >>",
                    quote(outfile),
                ]
            )
        )
    return cmds


def _write_script(cmds, script):
    """Write commands into a bash script"""
    f = open(script, "w+")
    for cmd in cmds:
        print(cmd, file=f)
    f.close()


def syscall(cmd, verbose=False):
    """Manage error handling when making syscalls"""
    if verbose:
        print("Running command:", cmd, flush=True)
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        print(
            "The following command failed with exit code",
            error.returncode,
            file=sys.stderr,
        )
        print(cmd, file=sys.stderr)
        print("\nThe output was:\n", file=sys.stderr)
        print(error.output.decode(), file=sys.stderr)
        raise Error("Error running command:", cmd)
    if verbose:
        print(decode(output))


def run_cmd(cmds, verbose=False):
    """Write and excute script"""
    tmpdir = tempfile.mkdtemp(prefix="tmp.", dir=os.getcwd())
    original_dir = os.getcwd()
    os.chdir(tmpdir)
    script = "run_jobs.sh"
    _write_script(cmds, script)
    syscall("bash " + script, verbose=verbose)
    os.chdir(original_dir)
    shutil.rmtree(tmpdir)
