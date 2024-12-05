import pandas as pd
import numpy as np


def scRNAcomparison(scRNAs, RNAs):
    matches = [scRNAs.header + "," + RNAs.header]
    for rna1 in scRNAs.sequences:
        for rna2 in RNAs.sequences:
            if rna1.gene == rna2.gene:
                matches.append([rna1.gene]+rna1.values+[rna2.gene]+rna2.values)
    pd.DataFrame(np.array(matches)).to_csv(
        "results/data/rnaSeqComparison___" + scRNAs.name + "___" + RNAs.name + ".csv",
        index=False, header=False)
