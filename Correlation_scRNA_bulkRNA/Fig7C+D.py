from libs.structure.RNASequencing import RNASequencing
from libs.structure.SingleCellRNASequencing import SingleCellRNASequencing
from libs.evaluation.singleCellRNA_brockenAxis import single_cell_RNA_bA
from libs.evaluation.singleCellRNA import single_cell_RNA

single_cell_rna_seq = SingleCellRNASequencing("raw_data/singleCellRNASeq/scRNA-seq_DEGs_WT_RARa-KO_RARa-TG.csv")

rnas = {"RARa-KO blank vs WT blank" : RNASequencing("raw_data/RNA_seq/RARa-KO blank vs WT blank.csv"),
        "BATF-KO blank vs WT blank" : RNASequencing("raw_data/RNA_seq/BATF-KO blank vs WT blank.csv")}

single_cell_RNA_bA(single_cell_rna_seq, rnas["RARa-KO blank vs WT blank"], "KO_Log2FoldChange")
single_cell_RNA(single_cell_rna_seq, rnas["BATF-KO blank vs WT blank"], "TG_Log2FoldChange")
