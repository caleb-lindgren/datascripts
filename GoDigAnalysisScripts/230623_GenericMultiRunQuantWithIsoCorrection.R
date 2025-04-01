library(rstudioapi)
library(ggplot2)
library(svMisc)
library(tidyr)
library(stringr)
library(dplyr)

#### 1. Parameters ####

# Isotopic correction table
isoTableLocation = "~/Research/TMT lot certificates of analysis/230414 XH364389/230414 XH364389.csv"

# Directory with GoDig output
wd = "D:/GoDig outputs/SRSIV012_SRSIII090_humanBrain-reAnalysis/"

# Files
filePattern = "^ed.*csv$" # Uses Regex

outputDirName = "250227_SRSIV012_humanBrainReanalysis"

prefix = ""

sumSN = 180

# Do you want to filter out all bad MS3 scans?
filterMS3wiseSumSN = T

# Don't enable this for site analysis unless the "gene symbols" specify the sites.
doGeneLevelAnalysis = T

#### 2. Save Parameters and Import ####

setwd(wd)
fileName = grep(filePattern, list.files(), value = T)
fileNames = paste(fileName, collapse = ", ")

outputDirExists = dir.exists(outputDirName)
dir.create(outputDirName)

destination = paste0(outputDirName, "\\", prefix)

params = c(wd, filePattern, fileNames, outputDirName, outputDirExists, prefix, destination, sumSN, filterMS3wiseSumSN, doGeneLevelAnalysis)
paramNames = c("Working directory", "Filename regex pattern", "Matching filename(s)", "Output directory", "Output directory already exists", "Output file prefix", "Whole output file path", "Sum SN Filter", "Filter for Sum SN MS3-wise", "Do Gene-Level Analysis")
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

#### 4. Isotope Correction Matrix ####

# Import
iso = read.csv(isoTableLocation)
names(iso)
names(iso) = c("MassTag", "ReporterIonMass", "Minus2C", "MinusCN", "MinusC", "NameMinusC", "MinusN", "NameMinusN", "M", "PlusN", "PlusC", "NamePlusC", "PlusNC", "Plus2C")

# Get rid of character columns
noNames = iso[,!grepl("Name", names(iso))]

# Change percent strings to real numbers
nn = noNames
for (i_col in 3:ncol(noNames)) {
  nn[,i_col] = gsub("\"", "", noNames[,i_col])
  nn[,i_col] = gsub("%", "", nn[,i_col])
  nn[,i_col] = as.numeric(nn[,i_col]) * 0.01
}

# Get the names back
nn[,grep("Name", names(iso), value = T)] = iso[,grep("Name", names(iso), value = T)]

# Start building the correction matrix

cm = matrix(rep(0, 18 * 18), nrow = 18, ncol = 18, 
            dimnames = list(paste0("channel_", 1:18), paste0("corr", 1:18)))

# Every channel (1-18, 126-135N) gets 100% contribution from itself
for (i_chan in 1:18)
  cm[i_chan,i_chan] = 1

# Odd channels 1-17 (126-134C) can have an additional 15N
channelInts = seq(1,17,2)
for (i_channel in channelInts) {
  cm[i_channel + 1, i_channel] = nn[i_channel, "PlusN"]
}

# All channels 1-16 (126-134N) can have an additional 13C
channelInts = 1:16
for (i_channel in channelInts) {
  cm[i_channel + 2, i_channel] = nn[i_channel, "PlusC"]
}

# Odd channels 1-15 (126-133C) can have an additional 13C AND 15N
channelInts = seq(1,15,2)
for (i_channel in channelInts) {
  cm[i_channel + 3,i_channel] = nn[i_channel, "PlusNC"]
}

# All channels 1-14 (126-133N) can have two additional 13C atoms
channelInts = 1:14
for (i_channel in channelInts) {
  cm[i_channel + 4,i_channel] = nn[i_channel, "Plus2C"]
}

# Even channels 2-18 (127N-135N) can lose a 15N
channelInts = seq(2,18,2)
for (i_channel in channelInts) {
  cm[i_channel - 1,i_channel] = nn[i_channel, "MinusN"]
}

# All channels 3-18 (127C-135N) can lose a 13C
channelInts = 3:18
for (i_channel in channelInts) {
  cm[i_channel - 2,i_channel] = nn[i_channel, "MinusC"]
}

# Even channels 4-18 (127N-135N) can lose a 15N AND 13C
channelInts = seq(4,18,2)
for (i_channel in channelInts) {
  cm[i_channel - 3,i_channel] = nn[i_channel, "MinusCN"]
}

# All channels 5-18 (128C-135N) can lose two 13Cs
channelInts = 5:18
for (i_channel in channelInts) {
  cm[i_channel - 4,i_channel] = nn[i_channel, "Minus2C"]
}

#### 5. Correct isotopes ####

corrected = precData

for (i_row in 1:nrow(corrected)) {
  quants = corrected[i_row, grep("Quant", names(corrected))]
  corrected[i_row, grep("Quant", names(corrected))] = solve(cm, quants)
}

# Save 
write.csv(corrected, paste0(destination, "4_precursorQuantCorrected.csv"), row.names = F)

# Turn negative values into 0s

corrected2 = corrected

for (i_row in 1:nrow(corrected2)) {
  quants = corrected2[i_row, grep("Quant", names(corrected2))]
  quants[quants < 0] = 0
  corrected2[i_row, grep("Quant", names(corrected2))] = quants
}

write.csv(corrected2, paste0(destination, "5_precQuantNoNeg.csv"), row.names = F)

#### 6. Plot ####

corrected = corrected2

for (i_row in 1:nrow(corrected)) {
  progress(i_row, max.value = nrow(corrected))
  df = data.frame(Channel = 1:18, SN = t(corrected[i_row, grepl("Quant", names(corrected))]))
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

#### 7. Widen data ####

cor2 = corrected[,names(corrected) != "PrecRun" &
                   names(corrected) != "SumSN" &
                   names(corrected) != "FileIndex"]
names(cor2) = gsub("Quant", "Channel", names(cor2))
cor2$RunName = cor2$RunName
wide = pivot_wider(cor2, names_from = "RunName", values_from = grep("Channel", names(cor2)))

write.csv(wide, paste0(destination, "7_corrected_wideFormat_noNeg.csv"), row.names = F)

#### 8. Gene level analysis ####

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
            SumSN = sum(SumSN)
  )
write.csv(geneLong, paste0(destination, "8_corrected_geneQuant.csv"), row.names = F)

# Plot
for (i_row in 1:nrow(geneLong)) {
  progress(i_row, max.value = nrow(geneLong))
  df = data.frame(Channel = 1:18, SN = t(geneLong[i_row, grepl("Quant", names(geneLong))]))
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
          