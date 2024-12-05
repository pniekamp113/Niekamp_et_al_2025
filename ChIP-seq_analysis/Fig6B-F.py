from libs.structure.ChiPSequencing import ChiPSeq
from libs.evaluation.overlappingGens import overlapping_gens


BATF_CTRL = ChiPSeq(["raw_data/Chip_seq/BATF Ctrl.csv"], name="BATF ctrl", color="red")
BATF_RA = ChiPSeq(["raw_data/Chip_seq/BATF RA.csv"], name="BATF RA", color="orange")
RARa_CTRL = ChiPSeq(["raw_data/Chip_seq/RARa Ctrl.csv"], name=r"RAR$\alpha$ ctrl", no_latex_name="RARa ctrl", color="blue")
RARa_RA = ChiPSeq(["raw_data/Chip_seq/RARa RA.csv"], name=r"RAR$\alpha$ RA", no_latex_name="RARa_RA", color="green")

overlapping_gens(RARa_CTRL, RARa_RA)
overlapping_gens(BATF_CTRL, BATF_RA)
overlapping_gens(RARa_CTRL, BATF_CTRL)
overlapping_gens(RARa_RA, BATF_RA) 
