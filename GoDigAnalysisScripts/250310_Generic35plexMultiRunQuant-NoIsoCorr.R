# NOTE: NO ISOTOPE CORRECTION #

library(rstudioapi)
library(ggplot2)
library(svMisc)
library(tidyr)
library(stringr)
library(dplyr)

#### 1. Parameters ####

# Directory with GoDig output
wd = "D:/GoDig outputs/SRSIV018_DBIA-35plex-GoDig-doseResponse/"

# Protein normalization factor CSV (must be copy-pasted from Protein Quant in Core to Excel then saved as CSV)
pNormFile = "D:/Core outputs/SRSIV018 DBIA GoDig 35plex/250319_SRSIV018A-FT_NormFactors.csv"

# Normalize by protein?
pNorm = T

# Files
filePattern = "^ed.*IV018A.*Tyr.*csv$" # Uses Regex

outputDirName = "250328_IV018A-TyrKinaseCys_processed"

prefix = ""

sumSN = 350

# Rearrange channels to match core?
goDig2CoreChannelLayout = T

# Do you want to filter out all bad MS3 scans?
filterMS3wiseSumSN = T

# Don't enable this for site analysis unless the "gene symbols" specify the sites.
doGeneLevelAnalysis = T

# Incorporate custom "GeneSymbol" entries (e.g., site info?)
customGenes = F
customGeneFile = "D:/GoDig target lists/250310_MH_SRSIV014_532-phos_35plex_GoDigTargetsWithSiteInfo.csv"

#### 2. Save Parameters and Import ####

setwd(wd)
fileName = grep(filePattern, list.files(), value = T)
fileNames = paste(fileName, collapse = ", ")

outputDirExists = dir.exists(outputDirName)
dir.create(outputDirName)

destination = paste0(outputDirName, "\\", prefix)

params = c(wd, pNormFile, pNorm, filePattern, fileNames, outputDirName, outputDirExists, prefix, destination, sumSN, filterMS3wiseSumSN, doGeneLevelAnalysis, customGenes, customGeneFile, goDig2CoreChannelLayout)
paramNames = c("Working directory", "Protein Normalization Factor File", "Normalize by Total Protein", "Filename regex pattern", "Matching filename(s)", "Output directory", "Output directory already exists", "Output file prefix", "Whole output file path", "Sum SN Filter", "Filter for Sum SN MS3-wise", "Do Gene-Level Analysis", "Incorporate Custom Gene Names", "Custom Gene File", "Rearrange 35plex channels to match Core")
paramDF = data.frame(Parameter.Name = paramNames, Parameter = params)
write.csv(paramDF, paste0(destination, "1_params.csv"), row.names = F)

data = NULL
for (i_file in 1:length(fileName)) {
  myFile = fileName[i_file]
  myRun = str_sub(myFile, 1, str_length(myFile) - 4)
  myData = read.csv(myFile)
  if (nrow(myData) <= 0)
    next
  myData$FileIndex = i_file
  myData$Filename = myFile
  myData$RunName = myRun
  if (i_file == 1)
    data = myData
  else
    data = rbind(data, myData)
}


# Remove unnecessary columns
data1 = data[,c("FileIndex", "Filename", "RunName", "GeneSymbol", "Peptide", "z", grep("Quant", names(data), value = T), "Sum.SN", "Iso.Spec")]

# Add precursor column 
data1$Precursor = paste(data1$Peptide, data1$z, sep="_")

# Merge in custom genes
if (customGenes) {
  cg = read.csv(customGeneFile)
  data1 = merge(data1, cg)
  data1$GeneSymbol = data1$Site
}

# Save
write.csv(data1, paste0(destination, "2_allRuns_selectColumns.csv"), row.names = F)

#### 3. Combine MS3s for each precursor ####

data1$PrecRun = paste(data1$Precursor, data1$FileIndex, sep = "_")

# Get combined (summed) quantities per channel
quantCols = grep("Quant", names(data1), value = T)
precs = unique(data1$PrecRun)
quantList = list()

for (i_prec in 1:length(precs)) {
  progress(i_prec, max.value = length(precs))
  prec = precs[i_prec]
  sub = data1[data1$PrecRun == prec,]
  if (filterMS3wiseSumSN) {
    sub = sub[sub$Sum.SN >= sumSN,]
  }
  theseQuants = c()
  for (quantCol in quantCols) {
    theseQuants = c(theseQuants, sum(sub[,quantCol]))
  }
  names(theseQuants) = quantCols
  quantList[[prec]] = theseQuants
}

quants = data.frame(do.call(rbind, quantList))
quants$PrecRun = row.names(quants)

# Add in gene names
infoDF = unique(data1[,c("GeneSymbol", "Peptide", "z", "Precursor", "PrecRun", "FileIndex", "RunName")])
precData = merge(infoDF, quants)

# Add in Sum SN
precData$SumSN = rowSums(precData[,grep("Quant", names(precData))])

# Remove insufficient Sum SN precursors
precData = precData[precData$SumSN >= sumSN,]

# Save
write.csv(precData, paste0(destination, "3_precursorQuantUncorrected.csv"), row.names = F)

# Rearrange columns if desired
if (goDig2CoreChannelLayout) {
  gdQuants = 1:35
  middleChannels = rep(c(2, 18, 3, 19), 7) +
    rep(2 * 1:7, each = 4)
  coreQuants = c(1, 2, 3, 19, 
                 middleChannels,
                 18, 34, 35)
  newColNames = paste0("Quant.", coreQuants)
  names(precData)[grep("^Quant\\.", names(precData))] =
    newColNames
  newColNamesSorted = paste0("Quant.", 1:35)
  nonQuantNames = names(precData)[!grepl("^Quant\\.", names(precData))]
  precData = precData[c(nonQuantNames, newColNamesSorted)]
}

#### 4. Normalize by total protein ####

if (pNorm) {
  pn = read.csv(pNormFile)
  qCols = grep("Quant", names(precData), value = T)
  for (i_q in 1:length(qCols)) {
    precData[,qCols[i_q]] = 
      precData[,qCols[i_q]] * pn$NormFactor[i_q]
  }
}

#### 5. Plot ####

corrected = precData

for (i_row in 1:nrow(corrected)) {
  progress(i_row, max.value = nrow(corrected))
  df = data.frame(Channel = 1:35, SN = t(corrected[i_row, grepl("Quant", names(corrected))]))
  names(df) = c("Channel", "Corrected SN")
  df$Channel = as.factor(df$Channel)
  gene = corrected[i_row, "GeneSymbol"]
  precString = paste0(corrected[i_row, "Peptide"], " (z=", corrected[i_row, "z"], ")")
  myMax = ceiling(max(df$`Corrected SN` / 10)) * 10
  tickMarkInterval = round(myMax / 10, -1)
  if (tickMarkInterval <= 0)
    tickMarkInterval = round(myMax/10)
  newMax = tickMarkInterval * 11
  p = ggplot(df, aes(x = Channel, y = `Corrected SN`, fill = Channel)) +
    geom_bar(stat = "identity", color = 'black', width = 0.5) +
    theme_classic() +
    ggtitle(paste0(gene, " - ", precString)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'None',
          panel.grid.major.y = element_line(color = "#BBBBBB")) +
    scale_y_continuous(breaks = seq(0, newMax, tickMarkInterval))
  p
  filename = paste0(destination, "6_bar_", gene, "_", corrected$Precursor[i_row], "__", corrected$RunName[i_row], ".pdf")
  newFilename = gsub("M\\*", "oM", filename)
  ggsave(newFilename, p, width = 8, height = 3.5)
  newFilename2 = gsub("\\.pdf", ".png", newFilename)
  ggsave(newFilename2, p, width = 8, height = 3.5)
}

#### 6. Widen data ####

cor2 = corrected[,names(corrected) != "PrecRun" &
                   names(corrected) != "SumSN" &
                   names(corrected) != "FileIndex"]
names(cor2) = gsub("Quant", "Channel", names(cor2))
cor2$RunName = cor2$RunName
wide = pivot_wider(cor2, names_from = "RunName", values_from = grep("Channel", names(cor2)))

write.csv(wide, paste0(destination, "7_corrected_wideFormat_noNeg.csv"), row.names = F)

#### 7. Gene level analysis ####

if (doGeneLevelAnalysis) {

# Analyze & save
geneLong = corrected %>% group_by(GeneSymbol) %>%
  summarize(Runs = length(unique(RunName)),
            Peptides = length(unique(Peptide)),
            Precursors = length(unique(Precursor)),
            Quant.1 = sum(Quant.1),
            Quant.2 = sum(Quant.2),
            Quant.3 = sum(Quant.3),
            Quant.4 = sum(Quant.4),
            Quant.5 = sum(Quant.5),
            Quant.6 = sum(Quant.6),
            Quant.7 = sum(Quant.7),
            Quant.8 = sum(Quant.8),
            Quant.9 = sum(Quant.9),
            Quant.10 = sum(Quant.10),
            Quant.11 = sum(Quant.11),
            Quant.12 = sum(Quant.12),
            Quant.13 = sum(Quant.13),
            Quant.14 = sum(Quant.14),
            Quant.15 = sum(Quant.15),
            Quant.16 = sum(Quant.16),
            Quant.17 = sum(Quant.17),
            Quant.18 = sum(Quant.18),
            Quant.19 = sum(Quant.19),
            Quant.20 = sum(Quant.20),
            Quant.21 = sum(Quant.21),
            Quant.22 = sum(Quant.22),
            Quant.23 = sum(Quant.23),
            Quant.24 = sum(Quant.24),
            Quant.25 = sum(Quant.25),
            Quant.26 = sum(Quant.26),
            Quant.27 = sum(Quant.27),
            Quant.28 = sum(Quant.28),
            Quant.29 = sum(Quant.29),
            Quant.30 = sum(Quant.30),
            Quant.31 = sum(Quant.31),
            Quant.32 = sum(Quant.32),
            Quant.33 = sum(Quant.33),
            Quant.34 = sum(Quant.34),
            Quant.35 = sum(Quant.35),
            SumSN = sum(SumSN)
  )
write.csv(geneLong, paste0(destination, "8_corrected_geneQuant.csv"), row.names = F)

# Plot
for (i_row in 1:nrow(geneLong)) {
  progress(i_row, max.value = nrow(geneLong))
  df = data.frame(Channel = 1:35, SN = t(geneLong[i_row, grepl("Quant", names(geneLong))]))
  names(df) = c("Channel", "Relative Abundance")
  df$Channel = as.factor(df$Channel)
  gene = geneLong[i_row, "GeneSymbol"]
  myMax = ceiling(max(df$`Relative Abundance` / 10)) * 10
  tickMarkInterval = round(myMax / 10, -1)
  if (tickMarkInterval <= 0)
    tickMarkInterval = round(myMax/10)
  newMax = tickMarkInterval * 11
  p = ggplot(df, aes(x = Channel, y = `Relative Abundance`, fill = Channel)) +
    geom_bar(stat = "identity", color = 'black', width = 0.5) +
    theme_classic() +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'None',
          panel.grid.major.y = element_line(color = "#BBBBBB")) +
    scale_y_continuous(breaks = seq(0, newMax, tickMarkInterval))
  p
  filename = paste0(destination, "9_bar_", gene, ".pdf")
  ggsave(filename, p, width = 8, height = 3.5)
  newFilename2 = gsub("\\.pdf", ".png", filename)
  ggsave(newFilename2, p, width = 8, height = 3.5)
}

}
          