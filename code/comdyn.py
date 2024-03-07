#!/usr/bin/env python

import argparse
import os
import re
import sys

import pandas as pd
from Bio.SeqUtils import seq1  # AA to oneletter notation
from pyfaidx import Fasta
import datetime

def parse_fasta_alignment(fasta_file="out.fasta", verbose=False, flip=False):
    """Parses the alignment in FASTA file

    Args:
        fasta_file (str, required): Fasta file. Defaults to "out.fasta".
        verbose (bool, optional): Complain. Defaults to False.
        flip (bool, optional): Flip alignment query/target. Defaults to False.

    Returns:
        dict: coordinate translation
    """
    # Load fasta file
    sequences = Fasta(fasta_file)

    # Grab keys:
    keys = list(sequences.keys())

    if flip:
        """Flip Target and Query order:"""
        keys.reverse()

    translation_dict = {}
    t_index = 0
    q_index = 0

    for (
        query,
        target,
    ) in zip(str(sequences[keys[0]]), str(sequences[keys[1]])):
        # Should be impossible state:
        if target == "-" and query == "-":
            print(f"impossible state at {q_index}->{t_index}", file=sys.stderr)
            continue

        # Match/Mismatch
        if query != "-" and target != "-":
            q_index += 1
            t_index += 1
            if verbose:
                print(f"{q_index} {query}-{target} {t_index}")
            translation_dict[q_index] = t_index
            continue

        # Insertion
        if target == "-":
            q_index += 1
            if verbose:
                print(f"{q_index} {query}-{target} {None}")
            translation_dict[q_index] = None
            continue

        # Deletion
        if query == "-":
            t_index += 1
            if verbose:
                print(f"{None} {query}-{target} {t_index}")
            # Don't care about this option so we don't store it in a dict.
            continue
    return translation_dict


def parse_atom_block(itp=None, verbose=False):
    """Parse atom block

    Args:
        itp (str, required): itp file. Defaults to None.
        verbose (bool, optional): complain. Defaults to False.

    Returns:
        DataFrame: with columns: "id", "type", "resnr", "residu", "atom", "cgnr", "charge".
    """
    fh = open(itp, "r")
    atom_section = False
    line_nr = 1
    atom_block = []
    for line in fh:
        if line.startswith("[ atoms ]"):
            if verbose:
                print(f"FIRST ATOM section linenr: {line_nr}")
            atom_section = True

        if atom_section:
            line = line.strip()
            fields = re.split(r" +", line)
            if fields[0].startswith("#"):
                pass
            elif len(line) == 0:
                if verbose:
                    print(f"LAST ATOM section linenr : {line_nr-1}")
                atom_section = False
            elif len(fields) > 6:
                if fields[4] == "BB":
                    fields[0] = int(fields[0])
                    fields[2] = int(fields[2])
                    fields[5] = int(fields[5])
                    fields[6] = float(fields[6])
                    atom_block.append(fields[0:7])
        line_nr += 1

    df = pd.DataFrame(atom_block)
    df.columns = ["id", "type", "resnr", "residu", "atom", "cgnr", "charge"]

    return df


def conserved_rubberbands(
    data_file1, data_file2, homology_dict, seq1, seq2, threshold=0.1, verbose=False
):
    """
    Get the rubber bands to keep.

    Args:
        data_file1 (DataFrame): Dataframe obtained from `parse_itp()`
        data_file2 (DataFrame): Dataframe obtained from `parse_itp()`
        homology_dict (dict): Coordinate translation dict obtained from `parse_fasta_alignment()`
        seq1 (str): Amino acid sequence of file1
        seq2 (str): Amino acid sequence of file2
        threshold (float, optional): Threshold upon rubber bands are not conserved (in nm). Defaults to 0.1.
        verbose (bool, optional): complain. Defaults to False.

    Returns:
        list: rubber bands to conserve
    """
    to_keep = {}
    to_skip = {}

    reasons = {}
    reasons["RB where AA not other file"] = 0
    reasons["RB in both files. Size diff < threshold"] = 0
    reasons["RB in both files. Size diff >= threshold"] = 0
    reasons["RB not in other file"] = 0

    for res_tuple, bbb1, bbb2, d1 in zip(
        data_file1["res_tuple"],
        data_file1["bbb1"],
        data_file1["bbb2"],
        data_file1["distance"],
    ):
        # Residues which span a rubber band. BBB are linked to two residues.
        resnr1_file1 = res_tuple[0]
        resnr2_file1 = res_tuple[1]

        # If a location is present in other protein:
        resnr1_file2 = None
        if resnr1_file1 in homology_dict:
            resnr1_file2 = homology_dict[resnr1_file1]

        # If b location is present in other protein:
        resnr2_file2 = None
        if resnr2_file1 in homology_dict:
            resnr2_file2 = homology_dict[resnr2_file1]

        ### Lookup AA residues:
        aa1_file1 = seq1[resnr1_file1-1]
        aa2_file1 = seq1[resnr2_file1-1]
        aa1_file2 = None
        aa2_file2 = None

        if resnr1_file2:
            aa1_file2 = seq2[resnr1_file2 - 1]
        if resnr2_file2:
            aa2_file2 = seq2[resnr2_file2 - 1]

        info = f"bbb:{bbb1},{bbb2}: {resnr1_file1}({aa1_file1})-{resnr2_file1}({aa2_file1}) == {resnr1_file2}({aa1_file2})-{resnr2_file2}({aa2_file2})"
        if not (resnr1_file2 or resnr2_file2):
            # FIXME residue from seq1 not used:
            reason = f"SKIP: {info}. RB where AA not other file"
            reasons["RB where AA not other file"] += 1
            to_skip[(bbb1, bbb2)] = reason
            if verbose:
                print(reason)
        else:
            # Filter condition:
            c = data_file2["res_tuple"] == (resnr1_file2, resnr2_file2)
            found_distances = data_file2["distance"][c]

            if len(found_distances) == 1:
                d2 = float(found_distances.iloc[0])
                diff = abs(d2 - d1)

                # So, rubberband is in both files, add tuple:
                if diff < threshold:
                    # Less than threshold, so add to keep:
                    reason = f"KEEP: {info}. d:{d2:.2f}-{d1:.2f}={diff:.2f}. diff < {threshold}"
                    reasons["RB in both files. Size diff < threshold"] += 1
                    if verbose:
                        print(reason)
                    to_keep[(bbb1, bbb2)] = reason
                else:
                    # Not less than threshold, so add to skip:
                    reason = f"SKIP: {info} d:{d2:.2f}-{d1:.2f}={diff:.2f}. diff >= {threshold}"
                    reasons["RB in both files. Size diff >= threshold"] += 1
                    if verbose:
                        print(reason)
                    to_skip[(bbb1, bbb2)] = reason
            else:
                if len(found_distances):
                    print(len(found_distances))
                # Rubberband not found in other file. Don't know how to compare; so keep:
                reasons["RB not in other file"] += 1
                reason = f"SKIP: {info}. Rubber band not found in other file"
                if verbose:
                    print(reason)
                to_skip[(bbb1, bbb2)] = reason

    print(
        f"To keep: {len(set(to_keep))}/{len(data_file1)} {len(set(to_keep))/len(data_file1)*100:.1f}%"
    )
    print(
        f"To skip: {len(set(to_skip))}/{len(data_file1)} {len(set(to_skip))/len(data_file1)*100:.1f}%"
    )
    print()
    for key in reasons:
        print(key, ":", reasons[key])

    return to_keep, to_skip


def parse_rubber_band_block(itp=None):
    """Given an ITP file, parses the rubber/elastic band block.

    Args:
        itp (str, required): itp file. Defaults to None.

    Returns:
        DataFrame: Full rubber band block. With fields: bbb1, bbb2, chain, distance.
    """
    fh = open(itp, "r")
    rubber_section = False
    line_nr = 0
    rubber_bands = []
    for line in fh:
        line_nr += 1
        if line.startswith("#"):
            continue

        if line.startswith("; Rubber band"):
            print(f"Start rubber band section: {line_nr}")
            rubber_section = True

        if rubber_section:
            line = line.strip()
            line = re.split(r" +", line)
            if line[0].startswith(";"):
                pass
            elif len(line) == 5:
                rubber_bands.append(
                    [int(line[0]), int(line[1]), int(line[2]), float(line[3])]
                )
            else:
                print(f"End rubber band section: {line_nr-1}")
                rubber_section = False
    rubber_bands = pd.DataFrame(rubber_bands)
    rubber_bands.columns = ["bbb1", "bbb2", "chain", "distance"]
    return rubber_bands


def parse_itp(itp_file, verbose=False):
    """Wrapper fucntion to parse the full ITP file

    Args:
        itp_file (str): itpfile
        verbose (bool, optional): complain. Defaults to False.

    Returns:
        DataFrame: Dataframe with fields: "bbb1" "bbb2"	"chain"	"distance" "resnr1" "resnr2" "residu1" "residu2"
        str: amino acid sequence
    """
    atom_block = parse_atom_block(itp=itp_file, verbose=verbose)
    aa_seq = seq1("".join(atom_block["residu"].to_list()))

    cg2resnr = dict(zip(atom_block["cgnr"], atom_block["resnr"]))
    cg2resname = dict(zip(atom_block["cgnr"], atom_block["residu"]))

    combined_data = parse_rubber_band_block(itp=itp_file)

    combined_data["resnr1"] = combined_data["bbb1"].apply(lambda x: cg2resnr[x])
    combined_data["resnr2"] = combined_data["bbb2"].apply(lambda x: cg2resnr[x])
    combined_data["residu1"] = combined_data["bbb1"].apply(lambda x: cg2resname[x])
    combined_data["residu2"] = combined_data["bbb2"].apply(lambda x: cg2resname[x])

    combined_data["res_tuple"] = combined_data.apply(
        lambda x: (x["resnr1"], x["resnr2"]), axis=1
    )

    # Add prefix if PDB does not start at the AA number1:
    if min(combined_data["resnr1"]) != 1:
        aa_seq = "X" * (min(combined_data["resnr1"]) - 1) + aa_seq

    return combined_data, aa_seq

def uniquify(path):
    """
    Creates a unique filename by appending a number to the base filename if it already exists.
    Args:  
        path (str, required): The original filename.
    Return: 
        str: A unique filename.
    """
    filename, extension = os.path.splitext(path)
    counter = 1

    # Check if the filename already contains a number in parentheses
    match = re.search(r"(\d+)$", filename)
    if match:
        existing_number = int(match.group(1))
        filename = filename[:match.start()]  # Remove the existing number

        # Increment the counter to start from the next number
        counter = existing_number + 1

    while os.path.exists(path):
        path = f"{filename}{extension}.{counter}"
        counter += 1

    return path

def write_itp(itp_in=None, itp_out="testje.itp", to_skip=None, verbose=False):
    """Write filtered itp file

    Args:
        itp_in (DataFrame, required): Input itp dataframe.
        itp_out (str, optional): Output itp file. Defaults to "testje.itp".
        to_skip (set of tuples, required): Rubberbands to remove.
        verbose (bool, optional): Complain. Defaults to False.
    """

    if os.path.isfile(itp_out):
        new_file = uniquify(itp_out)
        print(f"Making backup file of {itp_out} to {new_file}")
        
        os.rename(itp_out, new_file)


    print(f'\nWriting to "{itp_out}"...')
    written = 0
    removed = 0
    out_fh = open(itp_out, "w")
    fh = open(itp_in, "r")
    rubber_section = False
    line_nr = 0
    for line in fh:
        line_nr += 1
        if line.startswith("; Rubber band"):
            out_fh.write(line)
            written += 1

            if verbose:
                print(f"Start rubber band section: {line_nr}")
            rubber_section = True

            continue
        if rubber_section:
            line_strip = line.strip()
            fields = re.split(r" +", line_strip)
            # Comment lines. Don't do anything.
            if line and (line[0].startswith(";") or line[0].startswith("#")):
                out_fh.write(line)
                written += 1
                continue
            elif len(fields) == 5:
                location = (int(fields[0]), int(fields[1]))
                if location in to_skip:
                    reason = to_skip[location]
                    out_fh.write(f"; [{reason}] {line_strip}\n")
                    removed += 1
                    continue
                else:
                    out_fh.write(line)
                    written += 1
            else:
                if verbose:
                    print(f"End rubber band section: {line_nr-1}")
                rubber_section = False
        if not rubber_section:
            out_fh.write(line)
            written += 1

    print(f"WRITTEN {written} lines")
    print(f"SKIPPED {removed} lines")


def args_parser():
    """Build the argument parser

    Returns:
        ArgumentParser: argument parser
    """
    parser = argparse.ArgumentParser(
        prog="comdyn.py",
        description="Removes non-conserved rubber bonds from a GROMACS topology file based on homology of two protein sequences",
        epilog="Example: comdyn.py --itp_file1 data/4zw9/molecule_0.itp --itp_file2 data/5eqi/molecule_0.itp --nm 0.1",
    )

    parser.add_argument(
        "--itp_file1",
        help="First itp file.",
        type=str,
        dest="itp_file1",
        default="data/4zw9/molecule_0.itp",
    )

    parser.add_argument(
        "--itp_file2",
        help="second itp file.",
        type=str,
        dest="itp_file2"
    )
    parser.add_argument(
        "--itp_out",
        help="Output file. (Default: %(default)s)",
        type=str,
        dest="itp_out",
        default="output.itp",
    )
    parser.add_argument(
        "--nm",
        dest="nanometer",
        default=0.1,
        type=float,
        help="Nanometer cutoff to determine conserved rubber bands. (Default: %(default)s)",
    )
    parser.add_argument(
        "--gapopen",
        dest="gapopen",
        default=10.0,
        type=float,
        help="Gapopen penalty for NEEDLE alignment (Default: %(default)s)",
    )
    parser.add_argument(
        "--gapextend",
        dest="gapextend",
        default=0.5,
        type=float,
        help="Gapextend penalty for NEEDLE alignment (Default: %(default)s)",
    )
    parser.add_argument(
        "--verbose",
        dest="verbose",
        default=False,
        action="store_true",
        help="Complain. (Default: %(default)s)",
    )

    parser.add_argument(
        "--write",
        dest="write",
        default=False,
        action="store_true",
        help="Write new ITP files. (Default: %(default)s)",
    )
    parser.add_argument(
        "--inverse",
        dest="inverse",
        default=False,
        action="store_true",
        help="Keep the rubber bands that should be removed instead (Default: %(default)s)",
    )
    return parser


def print_args(args):
    """Print arguments

    Args:
        args (ArgumentParers): defined argument parser
    """
    print("-------------------------")
    print("Arguments:")
    print(f"itp file1 : {args.itp_file1}")
    print(f"itp file2 : {args.itp_file2}")
    print(f"Threshold : {args.nanometer}nm ({args.nanometer*10}A)")
    print("-------------------------")
    print("NEEDLE settings:")
    print(f"gapopen   : {args.gapopen}")
    print(f"gapextend : {args.gapextend}")
    print("output1   : out.fasta")
    print("output2   : out.needle")
    print("-------------------------")
    print(f"write itp : {args.write}")
    if args.write:
        print(f"itp out   : {args.itp_out}")
    print(f"verbose   : {args.verbose}")
    print("-------------------------")
    print()
    print()


def run_needle(gapopen=10.0, gapextend=0.5):
    """
    Performs an alignment using NEEDLE

    Args:
        gapopen (float, optional): gapopen parameter. Defaults to 10.
        gapextend (float, optional): gapextend parameter. Defaults to 0.5.
    """    

    try:
        if os.path.isfile("out.fasta"):
            print("Remove old alignment file 'out.fasta'")
            os.remove("out.fasta")
        
        print(
            f"needle -gapopen {gapopen} -gapextend {gapextend} -asequence a.fasta -bsequence b.fasta -outfile out.fasta -auto -aformat fasta"
        )
        os.system(
            f"needle -gapopen {gapopen} -gapextend {gapextend} -asequence a.fasta -bsequence b.fasta -outfile out.fasta -auto -aformat fasta"
        )
        os.system(
            f"needle -gapopen {gapopen} -gapextend {gapextend} -asequence a.fasta -bsequence b.fasta -outfile out.needle -auto"
        )
    except:
        sys.exit("FATAL: Some error in Needle alignment!")


def main():
    # Build argument parser
    parser = args_parser()
    args = parser.parse_args()
    print_args(args)

    # DO STUFF

    # Parse ITP files:
    print(f"Parsing \"{os.path.basename(args.itp_file1)}\"")
    itp_content_file1, aa_seq1 = parse_itp(
        itp_file=args.itp_file1,
        verbose=args.verbose,
    )

    print(f"\nParsing \"{os.path.basename(args.itp_file2)}\"")
    itp_content_file2, aa_seq2 = parse_itp(
        itp_file=args.itp_file2,
        verbose=args.verbose,
    )
    print()

    with open("a.fasta", "w") as f:
        f.write(">a\n")
        f.write(aa_seq1 + "\n")

    with open("b.fasta", "w") as f:
        f.write(">b\n")
        f.write(aa_seq2 + "\n")

    run_needle(gapopen=args.gapopen, gapextend=args.gapextend)
    print()
    homology_dict1 = parse_fasta_alignment(
        fasta_file="out.fasta", 
        flip=False, 
        verbose=args.verbose
    )

    # Calculate rubber bands to remove:
    if args.inverse:
        print("!!!   DOING THE OPPOSITE OF THIS !!!")
        print("!!! ITEMS TO KEEP are being REMOVED  !!!")

    # Conserved rubber bands:
    to_keep, to_skip = conserved_rubberbands(
        data_file1=itp_content_file1,  # From
        data_file2=itp_content_file2,  # To
        homology_dict=homology_dict1,
        seq1=aa_seq1,
        seq2=aa_seq2,
        threshold=args.nanometer,
        verbose=args.verbose,
    )

    print()

    if args.inverse:
        to_skip = to_keep

    if args.write:
        # Write:
        write_itp(
            itp_in=args.itp_file1,
            itp_out=args.itp_out,
            to_skip=to_skip,
            verbose=args.verbose,
        )
    else:
        print("\n!!!! NOT WRITING NEW ITP FILES !!!!")


if __name__ == "__main__":
    main()
