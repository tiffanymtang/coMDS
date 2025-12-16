data_dir <- here::here("data", "raw_data")
out_dir <- here::here("data")

############## cell trajectory ######################

data <- read.table(file.path(data_dir, "GSE98664_tpm_sailfish_mergedGTF_RamDA_mESC_differentiation_time_course.txt"), header = T)
rownames(data) <- data[, 1]
data <- data[, -1]

info <- colnames(data)
info <- gsub("RamDA_mESC_", "", info)
info <- substr(info, 1, 3)

data <- Seurat::CreateSeuratObject(counts = data, project = "TI", min.cells = 3, min.features = 200)
data <- Seurat::NormalizeData(data)

# identify highly variable genes
data <- Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = 500)
data <- Seurat::ScaleData(data, features = rownames(data))

# get the normalized data matrix
scale.data <- data@assays[["RNA"]]$scale.data
scale.data.var <- scale.data[which(data@assays[["RNA"]]@meta.data$var.features == rownames(scale.data)), ]

data <- t(scale.data.var)
save(data, info, file = file.path(out_dir, "traj.RData"))

############## cell cycle ######################

g2m <- read.table(file.path(data_dir, "E-MTAB-2805", "G2M_singlecells_counts.txt"), header = T)
s <- read.table(file.path(data_dir, "E-MTAB-2805", "S_singlecells_counts.txt"), header = T)
g1 <- read.table(file.path(data_dir, "E-MTAB-2805", "G1_singlecells_counts.txt"), header = T)

gene <- g1[, 1:4]
g1 <- g1[, -(1:4)]
s <- s[, -(1:4)]
g2m <- g2m[, -(1:4)]
combined <- cbind(g1, s, g2m)

# remove ERCCs
data <- combined[-which(substr(gene[, 1], 1, 4) == "ERCC"), ]
gene <- gene[-which(substr(gene[, 1], 1, 4) == "ERCC"), ]
info <- c(rep(c("G1", "S", "G2M"), each = 96))

# normalization and QC
data <- Seurat::CreateSeuratObject(counts = data, project = "cell_cycle", min.cells = 3, min.features = 20)
data <- Seurat::NormalizeData(data)

# identify highly variable genes
data <- Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- Seurat::ScaleData(data, features = rownames(data))

# get the normalized data matrix
scale.data <- data@assays[["RNA"]]$scale.data
scale.data.var <- scale.data[which(data@assays[["RNA"]]@meta.data$var.features == rownames(scale.data)), ]
data <- t(scale.data.var)

##### select relevant genes for cycling
cor.gene <- c()
for (i in 1:2000) {
  cor.gene[i] <- min(c(
    t.test(data[which(info == "G1"), i], data[-which(info == "G1"), i])$p.value,
    t.test(data[which(info == "G2M"), i], data[-which(info == "G2M"), i])$p.value,
    t.test(data[which(info == "S"), i], data[-which(info == "S"), i])$p.value
  ))
}
data <- data[, which(p.adjust(cor.gene, method = "BH") < 0.3)]
save(data, info, file = file.path(out_dir, "cycle.RData"))


############## four EQ ######################

data_4eq <- read.csv(file.path(data_dir, "4eq_pcadr_100d.csv"))
labels_4eq <- read.csv(file.path(data_dir, "4eq_labels.csv"))
data <- data_4eq[, -1]
info <- labels_4eq$x


set.seed(1)
subind1 <- sample(which(info == "b.cells"), 300)
subind2 <- sample(which(info == "cd14.monocytes"), 300)
subind3 <- sample(which(info == "naive.cytotoxic"), 300)
subind4 <- sample(which(info == "regulatory.t"), 300)

data <- data[c(subind1, subind2, subind3, subind4), ]
info <- info[c(subind1, subind2, subind3, subind4)]
save(data, info, file = file.path(out_dir, "4eq.RData"))

############## eight EQ ######################

data_8eq <- read.csv(file.path(data_dir, "8eq_pca.csv"))
label_8eq <- read.csv(file.path(data_dir, "8eq_labels.csv"))
data <- data_8eq[, -1]

info <- label_8eq$x

set.seed(1)

subind1 <- sample(which(info == "b.cells"), floor(499 / 3994 * 1600))
subind2 <- sample(which(info == "cd14.monocytes"), floor(600 / 3994 * 1600))
subind3 <- sample(which(info == "cd4.t.helper"), floor(400 / 3994 * 1600))
subind4 <- sample(which(info == "cd56.nk"), floor(600 / 3994 * 1600))
subind5 <- sample(which(info == "memory.t"), floor(500 / 3994 * 1600))
subind6 <- sample(which(info == "naive.cytotoxic"), floor(398 / 3994 * 1600))
subind7 <- sample(which(info == "naive.t"), floor(499 / 3994 * 1600))
subind8 <- sample(which(info == "regulatory.t"), floor(498 / 3994 * 1600))

data <- data[c(subind1, subind2, subind3, subind4, subind5, subind6, subind7, subind8), ]
info <- info[c(subind1, subind2, subind3, subind4, subind5, subind6, subind7, subind8)]

save(data, info, file = file.path(out_dir, "8eq.RData"))

############## HIV ######################

library(reticulate)
np <- import("numpy")

hiv_labels <- np$load(file.path(data_dir, "hiv_label.npy"), allow_pickle = TRUE)
hiv_data <- np$load(file.path(data_dir, "hiv_70.npy"))

set.seed(1)

subind1 <- sample(which(hiv_labels == "B cell"), 200)
subind2 <- sample(which(hiv_labels == "CTLs"), 200)
subind3 <- sample(which(hiv_labels == "DCs"), 200)
subind4 <- sample(which(hiv_labels == "monocyte"), 200)
subind5 <- sample(which(hiv_labels == "NK cell"), 200)
subind6 <- sample(which(hiv_labels == "plasmablast"), 200)
subind7 <- sample(which(hiv_labels == "T cell"), 200)

data <- hiv_data[c(subind1, subind2, subind3, subind4, subind5, subind6, subind7), ]
info <- hiv_labels[c(subind1, subind2, subind3, subind4, subind5, subind6, subind7)]
save(data, info, file = file.path(out_dir, "hiv.RData"))

############## star ######################

data_orig <- read.csv(file.path(data_dir, "star_classification.csv"), header = T)
data <- data_orig |>
  dplyr::select(alpha, delta, u, g, r, i, z, redshift, class) |>
  dplyr::filter(
    u != -9999,
    g != -9999,
    r != -9999,
    i != -9999,
    z != -9999
  )
set.seed(1)
data <- data[
  c(
    sample(which(data$class == "GALAXY"), 600),
    sample(which(data$class == "QSO"), 200),
    sample(which(data$class == "STAR"), 200)
  ),
]

info <- data$class
data <- data |>
  dplyr::select(-class) |>
  scale()

data <- data[, -(1:2)]
save(data, info, file = file.path(out_dir, "star_rs.RData"))

############## wholesale ######################

data_orig <- data.table::fread(
  file.path(file.path(data_dir, "Wholesale customers data.csv")),
  data.table = FALSE
)
data_orig[which(data_orig$Channel == 2), ]$Channel <- "Retail"
data_orig[which(data_orig$Channel == 1), ]$Channel <- "Hotel/Restaurant/Cafe"
data_orig[which(data_orig$Region == 1), ]$Region <- "Lisbon"
data_orig[which(data_orig$Region == 2), ]$Region <- "Oporto"
data_orig[which(data_orig$Region == 3), ]$Region <- "Other Region"

data <- data_orig |>
  dplyr::select(-Channel, -Region) |>
  dplyr::mutate(
    dplyr::across(
      tidyselect::everything(),
      ~ log(.x + 1)
    )
  ) |>
  scale()
info <- data_orig$Channel
save(data, info, file = file.path(out_dir, "wholesale.RData"))

############## wheat ######################

raw_data <- read.table(file.path(data_dir, "seeds_dataset.txt"), dec = "\t")
raw_data[which(raw_data$V8 == 1), ]$V8 <- "Kama"
raw_data[which(raw_data$V8 == 2), ]$V8 <- "Rosa"
raw_data[which(raw_data$V8 == 3), ]$V8 <- "Canadian"
data.nolabel <- raw_data[, -c(8)]
data.nolabel$V1 <- as.numeric(data.nolabel$V1)
data.nolabel$V2 <- as.numeric(data.nolabel$V2)
data.nolabel$V3 <- as.numeric(data.nolabel$V3)
data.nolabel$V4 <- as.numeric(data.nolabel$V4)
data.nolabel$V5 <- as.numeric(data.nolabel$V5)
data.nolabel$V6 <- as.numeric(data.nolabel$V6)
data.nolabel$V7 <- as.numeric(data.nolabel$V7)
data <- scale(data.nolabel)
info <- raw_data$V8
save(data, info, file = file.path(out_dir, "wheat.RData"))

############## olive ######################

data <- read.csv(file.path(data_dir, "olive_raw.csv"), header = T)
info <- data$X
data.nolabel <- data[, c(4:11)]
data <- scale(data.nolabel)
save(data, info, file = file.path(out_dir, "olive.RData"))
