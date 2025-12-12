# This script uses the NCRF output (GAA_summary_barcode95.summary) to generate a histogram of the repeat lengths.

#### input parameters ####

ID <- "barcode95" # the sample identifier
Locus <- "FGF14" # the gene symbol
Motif <- "GAA" # the repeat motif
Cutoff = 600 # has to be adjusted to separate the long and short allele

#### custom function to round repeat number ####

round_half_up <- function(x, digits = 0) {
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5
  z <- trunc(z)
  z / 10^digits * posneg
}

#### build histogram ####

  path_to_data <- paste0("C:/Users/Fienemann/Documents/PD_Panel_Paper/NCRF/", Motif, "_summary_", ID, ".summary")
  Dataframe <- read.table(path_to_data)
  colnames(Dataframe)<- c("#line","motif","seq", "start", "end", "strand", "seqLen", "querybp", "mRatio", "m", "mm", "i", "d")
  
  hist_title <- paste0(ID,", ", Locus, " repeat expansion")
  hist_unfiltered <- gghistogram(Dataframe$querybp, title = hist_title, xlab = "repeat length (bp)", ylab = "number of reads", fill = "#98D1FF",  bins = 30) +
    geom_vline(xintercept = Cutoff, linetype = "longdash")
  
  print(hist_unfiltered)
  
  min_bp <- min(Cutoff)
  max_bp <- max(Dataframe$querybp+1)
  Dataframe <- filter(Dataframe, querybp > min_bp)
  Dataframe <- filter(Dataframe, querybp < max_bp)
  
  median_rl <- median(Dataframe$querybp)
  print(median_rl)
  
  hist_title <- paste0(min_bp, " < repeat length < ", max_bp)
  hist_filtered <- gghistogram(Dataframe$querybp, title = hist_title, xlab = "repeat length (bp)", ylab = "number of reads", fill = "#98D1FF",  bins = 30) +
    
    annotate("text", x=850, y=4.5, label=paste0("median: ", median_rl), angle=0, size = 4) +
    annotate("text", x=850, y=4, label=paste0("motif:", Motif), angle=0, size = 4) +
    annotate("text", x=850, y=3.5, label=paste0("median RN: ", round_half_up(median_rl/nchar(Motif))), angle=0, size = 4)
  
  print(hist_filtered)

#### export histogram ####

output_name = paste0(ID, "_", Locus, "_histogram.tiff")
tiff(output_name, units="cm", width=30, height=10, res=300)
plot_grid(hist_unfiltered, hist_filtered, nrow = 1, labels = c('A', 'B'))
dev.off()