from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy


def overlapping_gens(ChiP1, ChiP2):
    overlapps = []
    nonMatchingPeaks1 = []
    nonMatchingPeaks2 = []

    peaks1 = deepcopy(ChiP1.peaks)
    peaks2 = deepcopy(ChiP2.peaks)
    for i in range(len(peaks1)):
        peak1 = peaks1[i]
        save = True
        for j in range(len(peaks2)):
            peak2 = peaks2[j]
            if peak1 == peak2:
                save = False
                del peaks2[j]
                break
        if save:
            nonMatchingPeaks1.append([peak1.gene, peak1.fold_enrichment, peak1.peak_start, peak1.peak_end, peak2.fold_enrichment, peak2.peak_start, peak2.peak_end])

    peaks1 = deepcopy(ChiP1.peaks)
    peaks2 = deepcopy(ChiP2.peaks)
    for i in range(len(peaks2)):
        peak1 = peaks2[i]
        save = True
        for j in range(len(peaks1)):
            peak2 = peaks1[j]
            if peak1 == peak2:
                save = False
                del peaks1[j]
                break
        if save:
            nonMatchingPeaks2.append([peak1.gene, peak1.fold_enrichment, peak1.peak_start, peak1.peak_end, peak2.fold_enrichment, peak2.peak_start, peak2.peak_end])

    peaks1 = deepcopy(ChiP1.peaks)
    peaks2 = deepcopy(ChiP2.peaks)
    for i in range(len(peaks1)):
        peak1 = peaks1[i]
        for j in range(len(peaks2)):
            peak2 = peaks2[j]
            if peak1 == peak2:
                overlapps.append([peak1.gene, peak1.fold_enrichment, peak1.peak_start, peak1.peak_end, peak2.fold_enrichment, peak2.peak_start, peak2.peak_end])
                del peaks2[j]
                break

    number_overlapps = len(overlapps)
    venn2(subsets={'10': ChiP1.number_gens-number_overlapps, '01': ChiP2.number_gens-number_overlapps, '11': number_overlapps},
          set_labels=(ChiP1.name, ChiP2.name), set_colors=(ChiP1.color, ChiP2.color))
    plt.tight_layout()
    plt.savefig("results/images/overlappingGens___"+str(ChiP1.MAXOVERLAPP)+"___"+ChiP1.no_latex_name+"___"+ChiP2.no_latex_name+".png")
    plt.close()

    pd.DataFrame(overlapps).to_csv("results/data/overlappingGens___"+str(ChiP1.MAXOVERLAPP)+"___"+ChiP1.no_latex_name+"___"+ChiP2.no_latex_name+ ".csv", index=False, header=False)
    pd.DataFrame(nonMatchingPeaks1).to_csv(
        "results/data/overlappingGens___"+"only"+str(ChiP1.MAXOVERLAPP)+"___"+ ChiP1.no_latex_name+"___" + ChiP1.no_latex_name + "___" + ChiP2.no_latex_name + ".csv", index=False,
        header=False)
    pd.DataFrame(nonMatchingPeaks2).to_csv(
        "results/data/overlappingGens___"+"only"+str(ChiP1.MAXOVERLAPP)+"___"+ ChiP2.no_latex_name+"___" + ChiP1.no_latex_name + "___" + ChiP2.no_latex_name + ".csv", index=False,
        header=False)
