#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd

def main():
    df = []
    for protein in ["4zw9", "5eqi"]:
        for type in ["cg_crystal", "cg_equil", "cg_minim"]:
            for distance in ["TM1_TM7", "TM2_TM8", "TM5_TM11"]:
                for in_out in ["in", "out"]:
                    path = f"distances/{protein}/{type}_{distance}_{in_out}.com.xvg"
                    x, y = np.loadtxt(path, comments=["@", "#"], unpack=True)
                    description = f"{type}_{distance}_{in_out}"
                    res = {"protein" : protein,
                           "description" : description,
                           "distance"  : y}
                    df.append(res)
    df = pd.DataFrame(df)
    df.to_csv("initial_distances.tsv", sep="\t", index=False)
if __name__ == "__main__":
    main()
