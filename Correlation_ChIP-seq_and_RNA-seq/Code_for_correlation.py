import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import pandas as pd
import numpy as np
import matplotlib.colors as colors


def chip_vs_rna(chip, rna):
    fold_changes, fold_enrichment, color = [], [], []
    name = []
    for rna_seq in rna.sequences:
        for chip_gen in chip.peaks:
            if rna_seq.gene == chip_gen.gene:
                fold_changes.append(rna_seq.fold_change)
                fold_enrichment.append(chip_gen.fold_enrichment)
                color.append(rna_seq.p_value)
                name.append(rna_seq.gene)
    sc = plt.scatter(fold_changes, fold_enrichment, c=color, norm=colors.LogNorm(), s=20, cmap="coolwarm")
    plt.colorbar(sc).set_label("log p value")

    m, b, r_value, p_value, std_err = linregress(fold_changes, fold_enrichment)
    x = np.linspace(np.amin(fold_changes), np.amax(fold_changes), 10)
    plt.plot(x, m * x + b, color="black",
             label="slope: " + str(np.round(m, decimals=3)) + ", intersection:" + str(np.round(b, decimals=1))
                   + ", r_value: " + str(np.round(r_value, decimals=4)) + ", p_values: " + str(p_value))
    plt.legend(bbox_to_anchor=(1.3, 1.2))

    plt.xlabel(rna.name)
    plt.ylabel(chip.name)
    plt.subplots_adjust(top=0.8)
    plt.savefig("results/images/ChiPVsRNA___" + chip.no_latex_name + "___" + rna.name + ".png")
    plt.close()

    pd.DataFrame(np.array([name, fold_changes, fold_enrichment, color]).T).to_csv("results/data/ChiPVsRNA___" + chip.no_latex_name + "___" + rna.name + ".csv", index=False, header=False)
