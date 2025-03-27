#' @title plotbar_TopPerSample
#'
#' @description
#' This Function plots a bar with a top(i) is per sample instead of an overall overall top(i).
#' Idea is to have higher explained reads compared to the ps_plotbar which uses overall top abundance. The functions returns a GGPLOT object.
#' The top(i) is calculated in the taxrank, i.e if top = 5 and taxrank = "Family" the 5 most abundant families were taken. The taxfill is the color.The taxfill should
#' be equal or higher than the taxrank
#'
#' @keywords plotbar
#' @param ps phyloseq object to be used
#' @param top the amount of most abundant ranks to be added in the plot
#' @param relatief Is the graph y-axis with realtive or absolute abundance (TRUE or FALSE)
#' @param taxrank the rank for the top. Mostly Species is used (is default)
#' @param taxfill the color in the graph
#' @param x the x-axis names. If unspecified the ps object names will be used
#' @param legend.position: Where is the legend? (left, right, bottom, top, or none). Legends can be very large
#'
#' @param choosing "none" is equal to set legend=FALSE
#' @param GS Make a GenusSpecies rank in the ps? (TRUE or FALSE). Default to TRUE, makes a GenusSpecies rank
#' @param statistics: Do you want (relative) reads to be plotted (TRUE/FALSE)
#' @param angle angle for x-axis label (most likely 0 or 90°)
#' @param output: undocumented feature what to export graph or table the graph is made of
#' @export
#' @examples plotbar_TopPerSample(ps, taxfill="Family", taxrank="Genus", top=10)

#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 alpha
#' @importFrom stringr word
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq sample_sums


plotbar_TopPerSample <- function (ps, top = 10, relatief = TRUE, taxrank = "Species",
          taxfill = "Genus", output = "graph", x = "x_names", legend.position="right",
          statistics = FALSE, GS = TRUE, angle = 0) {

#  DEFAULT for Testing, can be removed
#  ps <- readRDS("c:/software/Rproject/MicrobiomeAdds/Data/tank_milk.Rds")
#  top <- 10
#  relatief <- FALSE
#  statistics <- TRUE
#  taxrank <- "Species"
#  taxfill <- "Genus"
#  output <- "graph"
#  x <- "x_names"
# legend.position="right"
#  GS <- TRUE
#  angle <- 90


# library(ggplot2)
#  library(microbiome)
#  library(microbiomeutilities)
#  library(tidyverse)
#  library(stringr)
#  library(data.table)

  ## first replace the TAG
  ## ps <- add_refseq(ps, tag = "ASV")
  ## extract sample table and add library size in reads
  samples <- data.frame(sample_data(ps)) %>%
    tibble::rownames_to_column(var = "x_naam") %>%
    mutate(x_name = x_naam) %>% tibble::column_to_rownames("x_naam")
  samples$LibSize = sample_sums(ps)

  ## mutate the sample name, or the row_name or the given x name value
  xvar <- sym(x)
  if (x != "x_names") {
    samples <- samples %>% mutate(x_name = !!xvar)
  }

  ##
  alle_samples <- unique(samples$x_name)
  print(paste("No. of samples:", length(alle_samples)))
  text_frame <- data.frame(x = character(0), label = character(0))

  ## add a GS level if not available
  TAX_tabel <- data.frame(tax_table(ps))
  if (GS) {
    TAX_tabel$GS_word <- word(TAX_tabel$Species, 1, sep = "\\(")
    TAX_tabel$GenusSpecies <- paste(TAX_tabel$Genus, TAX_tabel$GS_word)
    TAX_tabel <- TAX_tabel %>% rownames_to_column("rijnaam") %>%
      dplyr::select(rijnaam, Kingdom, Phylum, Class, Order,
                    Family, Genus, GenusSpecies, Species)
    start_col <- 10
  } else {
    TAX_tabel <- TAX_tabel %>% rownames_to_column("rijnaam")
    start_col <- 9
  }

  ## omzeilen van Otab <- data.frame(otu_table(ps)) --> erg traag! dit is (veel) sneller
  print("Extract OTU_tabel from Phyloseq")
  OTab <- otu_table(ps)
  OTab <- as.matrix(OTab@.Data)
  OTab <- data.frame(OTab)

  ## Connect TAX en de OTU tabellen uit de phyloseq
  OTU_tabel_T <- data.table::transpose(OTab, keep = "rijnaam")
  colnames(OTU_tabel_T) <- c("rijnaam", alle_samples)
  TAX_OTU <- left_join(TAX_tabel, OTU_tabel_T, by = "rijnaam") %>%
    mutate(sum = rowSums(across(where(is.numeric)))) %>%
    filter(sum >= 1) %>% dplyr::select(-sum)

  ## tax_glom for taxrank given
  glom_level <- taxrank
  TAX_OTU_A <- TAX_OTU
  if (glom_level == "Kingdom") {level_nr <- 2}
  if (glom_level == "Phylum")  {level_nr <- 3}
  if (glom_level == "Class")   {level_nr <- 4}
  if (glom_level == "Order")   {level_nr <- 5}
  if (glom_level == "Family")  {level_nr <- 6}
  if (glom_level == "Genus")   {level_nr <- 7}
  if (glom_level == "GenusSpecies") {level_nr <- 8}
  if (glom_level == "Species") {level_nr <- 9}

  ## level below taxrank then set to "<>"
  for (i in level_nr:9) {
    if (level_nr != i) {
      TAX_OTU_A[, i] <- "<>"
    }
  }

  TAX_OTU_sel <- TAX_OTU_A %>%
    mutate(KPCOFGgsS = paste0(Kingdom, "*", Phylum, "*", Class, "*", Order, "*", Family, "*",
                              Genus, "*", GenusSpecies, "*", Species)) %>%
    dplyr::select(ncol(TAX_OTU) + 1, 10:ncol(TAX_OTU)) %>%
    pivot_longer(names_to = "SampleName", values_to = "reads", cols = 2:(ncol(TAX_OTU) - 8)) %>%
    group_by(TAX_SN = paste(KPCOFGgsS, SampleName, sep = fixed("@"))) %>%
    summarize(reads = sum(reads)) %>%
    mutate(TAX = word(TAX_SN, 1, sep = fixed("@"))) %>%
    mutate(SampleName = word(TAX_SN, 2, sep = fixed("@"))) %>%
    dplyr::select(-TAX_SN)

  keer <- 1
  sample_names <- unique(TAX_OTU_sel$SampleName)
  tabel_stat <- data.frame(sample_name=character(0), reads.top=numeric(0), raeds.all=numeric(0), percentage=numeric(0))
  max_top <- 0

  ## for all samples: make a subset, arrange -reads, and select the given top. Stick togheter in uitvoer
  for (i in 1:length(sample_names)) {
    TAX_OTU_sel_sample <- subset(TAX_OTU_sel, TAX_OTU_sel$SampleName == sample_names[i])
    TAX_OTU_sel_sample <- TAX_OTU_sel_sample %>% arrange(-reads)
    TAX_OTU_sel_sample <- TAX_OTU_sel_sample[1:top, ]

    som_top <- sum(TAX_OTU_sel_sample$reads)
    if (max_top<som_top){max_top <- som_top} #find maximum sum top x reads

    x <- subset(samples, samples$x_name == TAX_OTU_sel_sample$SampleName[1])
    libSize = x$LibSize

    regel_stat <- cbind.data.frame(SampleName=TAX_OTU_sel_sample$SampleName[1], reads_in_top=som_top, allreads=libSize, percentage=round(som_top/libSize*100,1))
    tabel_stat <- rbind(tabel_stat, regel_stat)

    if (keer == 1) {
      uitvoer <- TAX_OTU_sel_sample

      keer <- 2
    } else {
      uitvoer <- rbind(uitvoer, TAX_OTU_sel_sample)
    }
  }
  uitvoer <- uitvoer %>%
    mutate(Kingdom = word(TAX, 1, sep = fixed("*"))) %>%
    mutate(Phylum =  word(TAX, 2, sep = fixed("*"))) %>%
    mutate(Class =   word(TAX, 3, sep = fixed("*"))) %>%
    mutate(Order =   word(TAX, 4, sep = fixed("*"))) %>%
    mutate(Family=   word(TAX, 5, sep = fixed("*"))) %>%
    mutate(Genus =   word(TAX, 6, sep = fixed("*"))) %>%
    mutate(GenusSpecies = word(TAX, 7, sep = fixed("*"))) %>%
    mutate(Species = word(TAX, 8, sep = fixed("*"))) %>%
    dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, GenusSpecies, Species, SampleName, reads)

  ## make wider
  uit_breed <- pivot_wider(uitvoer, names_from = "SampleName",
                           values_from = "reads")

  uitvoer2 <- left_join(uitvoer, tabel_stat, by="SampleName")
  uitvoer2$relatief <- relatief  ##relative or absolute reads?
  uitvoer3 <- uitvoer2 %>%
    mutate(y_value = ifelse(relatief == TRUE, reads/allreads*100, reads)) %>%
    mutate(stat_text = ifelse(relatief == TRUE, paste0(percentage,"%"), paste(reads_in_top,"/",allreads)))

  pl <- ggplot(uitvoer3, aes(x = SampleName, y = y_value, fill = fct_reorder(!!sym(taxfill), reads))) +
    theme_minimal() +
    geom_col(col = "black") +
    ylab("Abundance") +
    labs(fill=paste(sym(taxfill),"{fill color}"))+
    theme(axis.text.x = element_text(angle = angle, vjust = 0.5, hjust = 1))

  if (statistics == TRUE) {
    max_y = max(uitvoer3$y_value)
    if (max_y<=100){
        max_y<- 100
      } else {
        max_y <- max_top
        print(paste("Max y-value",max_y))
      }  ## if relatief then 100% scale set
    pl <- pl + coord_cartesian(ylim=c(0, max_y*1.2))
    pl <- pl +
      geom_text(data=uitvoer3, aes(label=stat_text, x=SampleName, y=max_y * 1.15, angle=90), size=3)
  } else {
    print("No statistics added")
  }

  ## The legend legends can be very large and preventing to see a graph by setting to "none"
  pl <- pl +
    theme(legend.position = legend.position)


  if (relatief==TRUE){
    pl <- pl + ylab("Relative Abundance (%)")
  } else {
      pl <- pl + ylab("Absolute Abundance (reads)")
  }

  ## output a plot or a table where the plot is made off
  if (output == "table") {
    pl <- uit_breed
  }
  return(pl)
}
