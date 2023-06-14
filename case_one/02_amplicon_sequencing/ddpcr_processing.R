## Most of the processing was done within the Biorad software. There are 3 files total, Section Plate 1 (Mut only).csv, Section Plate 2 (Mut only).csv and Final Analysis.csv. These are the two original plates and then Final Analysis was a table that the Biorad reps Chad and Camille gave me in which they tweaked the setting and output the file in a usable format from plate 1. For this reason I used Final Analysis for plate 1 and I had no choice but to use the original plate for plate 2:
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Digital_droplet_PCR")
datDDP1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Digital_droplet_PCR/Final Analysis.csv")
datDDP2 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Digital_droplet_PCR/Section Plate 2(Mut only).csv")
dim(datDDP1)
# [1] 71  6
dim(datDDP2)
# [1] 71 62
## Subset the ddPCR data from the second plate down so that both tables have the same columns:
datDDP2 = datDDP2[, is.element(colnames(datDDP2), colnames(datDDP1))]
dim(datDDP2)
# [1] 71  6
## There are blood and a NTC controls on each of the plates. I am only concerned with the sections as there was no mutation detected in each of these controls so take first 69 samples:
## Now just get the samples:
datDDP1Sub = datDDP1[1:69, ]
datDDP2Sub = datDDP2[1:69, ]
## Now rename the samples in the digital droplet samples:
datDDP1Sub$Sample = gsub("Section ", "", datDDP1Sub$Sample)
datDDP2Sub$Sample = gsub("Section ", "", datDDP2Sub$Sample)
## Convert the Sample to a numeric value:
datDDP1Sub$Sample = as.numeric(datDDP1Sub$Sample)
datDDP2Sub$Sample = as.numeric(datDDP2Sub$Sample)
## Order by section number:
datDDP1Sub = datDDP1Sub[order(datDDP1Sub$Sample), ]
datDDP2Sub = datDDP2Sub[order(datDDP2Sub$Sample), ]
## Now renumber these sequentially so that they match the renumbered gene expression data:
datDDP1Sub$Sample = c(1:69)
datDDP2Sub$Sample = c(1:69)
## Now add back the blood and NTC samples and write to disk:
datDDP1 = rbind(datDDP1Sub, datDDP1[70:71, ])
datDDP2 = rbind(datDDP2Sub, datDDP2[70:71, ])
write.table(datDDP1, file = "DDPCR_IDH1R132H_plate1_renumbered.csv", col.names = TRUE, row.names = FALSE, sep = ",")
write.table(datDDP2, file = "DDPCR_IDH1R132H_plate2_renumbered.csv", col.names = TRUE, row.names = FALSE, sep = ",")
