## In R combine the AmpSeq data from the first and second runs, using the IDH1 and ACCS from the first run that had higher coverage.
# {{{
## Read in the two data frames
amp1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/R_analysis_of_VarFreq/varFreq_first_batch_mutations_using_countreads.csv")
amp2 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Aligned_files_bwa/R_analysis_of_VarFreq/varFreq_second_batch_mutations_using_countreads.csv")
dim(amp1)
# [1] 70 73
dim(amp2)
# [1]  70 201
## First check that the sections are in the correct order
all.equal(amp1$sample, amp2$sample)
# [1] TRUE
## Now combine the data frames:
datCombi = data.frame(amp1, amp2[, 2:185])
## At this point I basically have the same table as I previously generated for the ampseq data but the colnames are different and the columns are in a slightly different order. To make replication of the previous analysis more straightforward I will rename the columns and reorder them.
## First write out the table with the readcounts etc present in it.
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Reanalysis_using_readcounts")
write.table(datCombi, file = "Combined_variant_frequencies_run1_and_run2_with_readCounts.csv", col.names = TRUE, row.names = FALSE, sep = ",")
## Now subset the data frame down to exclude the read counts.
datCombi = datCombi[, -grep("Q30.depth", colnames(datCombi))]
datCombi = datCombi[, -grep("Q30.reads.ref", colnames(datCombi))]
datCombi = datCombi[, -grep("Q30.reads.alt", colnames(datCombi))]
## I am still left with colnames that give the position in the amplicon whereas previously I used protein change as label. I will now relabel the columns appropriately.
colnames(datCombi) = gsub("IDH1.242.GtoA_", "IDH1.R132H.", colnames(datCombi))
colnames(datCombi) = gsub("TP53.137.TtoC_", "TP53.L145P.", colnames(datCombi))
colnames(datCombi) = gsub("ACCS.159.GtoA_", "ACCS.A197T.", colnames(datCombi))
colnames(datCombi) = gsub("GDPD2.214.GtoA_", "GDPD2.A104T.", colnames(datCombi))
colnames(datCombi) = gsub("FAM193B.223.CtoG_", "FAM193B.P556A.", colnames(datCombi))
colnames(datCombi) = gsub("CTNNA3.289.AtoG_", "CTNNA3.D364G.", colnames(datCombi))
colnames(datCombi) = gsub("TEC.120.TtoC_", "TEC.F92L.", colnames(datCombi))
colnames(datCombi) = gsub("GPR173.224.GtoA_", "GPR173.A41T.", colnames(datCombi))
colnames(datCombi) = gsub("IAPP.318.AtoG_", "IAPP.N55S.", colnames(datCombi))
colnames(datCombi) = gsub("ATRX.171.AtoDEL.3.GGA_", "ATRX.p.1459_1460del.", colnames(datCombi))
colnames(datCombi) = gsub("C2.152.GtoA_", "C2.R39Q.", colnames(datCombi))
colnames(datCombi) = gsub("CSAG1.101.GtoDEL.1.A_", "CSAG1.p.E64fs.", colnames(datCombi))
colnames(datCombi) = gsub("CTNND2.67.GtoDEL.3.GAA_", "CTNND2.p.814_815del.", colnames(datCombi))
colnames(datCombi) = gsub("DBC1.168.GtoA_", "DBC1.R622Q.", colnames(datCombi))
colnames(datCombi) = gsub("GPR126.277.GtoT_", "GPR126.C41F.", colnames(datCombi))
colnames(datCombi) = gsub("HIST1H4K.243.GtoA_", "HIST1H4K.G102D.", colnames(datCombi))
colnames(datCombi) = gsub("IL7R.261.GtoDEL.1.A_", "IL7R.p.K265fs.", colnames(datCombi))
colnames(datCombi) = gsub("MUC16.262.AtoT_", "MUC16.I1060F.", colnames(datCombi))
colnames(datCombi) = gsub("NTRK1.212.AtoG_", "NTRK1.E413G.", colnames(datCombi))
colnames(datCombi) = gsub("OR52N2.223.GtoA_", "OR52N2.R236H.", colnames(datCombi))
colnames(datCombi) = gsub("OR7C2.241.AtoDEL.1.T_", "OR7C2.p.F102fs.", colnames(datCombi))
colnames(datCombi) = gsub("PASD1.202.AtoG_", "PASD1.T679A.", colnames(datCombi))
colnames(datCombi) = gsub("PCLO.235.TtoC_", "PCLO.S506P.", colnames(datCombi))
colnames(datCombi) = gsub("PHF8.214.AtoC_", "PHF8.Q929P.", colnames(datCombi))
colnames(datCombi) = gsub("PLA2G4D.152.CtoT_", "PLA2G4D.A692V.", colnames(datCombi))
colnames(datCombi) = gsub("RECQL.212.AtoINS.1.A_", "RECQL.p.K40fs.", colnames(datCombi))
colnames(datCombi) = gsub("SLC20A2.97.GtoA_", "SLC20A2.R611H.", colnames(datCombi))
colnames(datCombi) = gsub("SUSD4.177.CtoT_", "SUSD4.Q37X.", colnames(datCombi))
colnames(datCombi) = gsub("TBCC.126.GtoT_", "TBCC.R12S.", colnames(datCombi))
colnames(datCombi) = gsub("WASL.245.TtoDEL.3.CCT_", "WASL.p.303_303del.", colnames(datCombi))
colnames(datCombi) = gsub("ZFYVE16.266.GtoA_", "ZFYVE16.D874N.", colnames(datCombi))
## Now change the other parts of the column names to match the original run.
colnames(datCombi) = gsub("sample", "section.no", colnames(datCombi))
colnames(datCombi) = gsub("freq", "var.freq", colnames(datCombi))
colnames(datCombi) = gsub("LCI", "L95CI", colnames(datCombi))
colnames(datCombi) = gsub("UCI", "U95CI", colnames(datCombi))
colnames(datCombi) = gsub("UCI", "U95CI", colnames(datCombi))
## It is clear that a number of the mutations were not detected in a single read in a single sample. These will be excluded from further analysis:
toexclude = c("ATRX", "CTNND2", "RECQL", "WASL")
datCombi = datCombi[, -grep("ATRX", colnames(datCombi))]
datCombi = datCombi[, -grep("CTNND2", colnames(datCombi))]
datCombi = datCombi[, -grep("RECQL", colnames(datCombi))]
datCombi = datCombi[, -grep("WASL", colnames(datCombi))]
dim(datCombi)
# [1]  70 141
## Write this table out to disk
write.table(datCombi, file = "Combined_variant_frequencies_run1_and_run2_without_readCounts_relabelled.csv", col.names = TRUE, row.names = FALSE, sep = ",")
# }}}
