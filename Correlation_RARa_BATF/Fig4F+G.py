from libs.structure.SingleCellRNASequencing import SingleCellRNASequencingGeneral
from libs.structure.RNASequencing import RNASequencingGeneral
from libs.evaluation.scRNAcomparison import scRNAcomparison

scRNAseq_RARaKO_vs_RARaTG = SingleCellRNASequencingGeneral("raw_data/singleCellRNASeq/scRNAseq_RARaKO_vs_RARaTG.csv")

Seo = RNASequencingGeneral("raw_data/singleCellRNASeq/Seo BATF overexpression vs control41590_2021_964_MOESM3_ESM.csv")
Tsao = RNASequencingGeneral("raw_data/singleCellRNASeq/Tsao 2021_BATF-KO gene regulation_sciimmunol.abi4919_table_s1.csv")

scRNAcomparison(scRNAseq_RARaKO_vs_RARaTG, Seo)
scRNAcomparison(scRNAseq_RARaKO_vs_RARaTG, Tsao)
