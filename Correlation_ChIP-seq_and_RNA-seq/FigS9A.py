from libs.structure.ChiPSequencing import ChiPSeq
from libs.structure.RNASequencing import RNASequencing
from libs.evaluation.ChipVsRNA import chip_vs_rna


rna = RNASequencing("raw_data/RNA_seq/WT RA vs Blank for chip seq comparison.csv")
RARa_CTRL = ChiPSeq(["raw_data/Chip_seq/RARa Ctrl.csv"], name=r"RAR$\alpha$ ctrl", no_latex_name="RARa ctrl", color="blue")

chip_vs_rna(RARa_CTRL, rna)
