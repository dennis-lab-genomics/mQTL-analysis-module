library(MatrixEQTL)
library(R.utils)

argv <- commandArgs(
  asValues = TRUE,
  defaults = list(
    data_dir = "",
    SNP_fname = "all_imputed_snps_matrixeQTL.txt",
    methylation_fname = "methylation_data_matrixeQTL.txt",
    cov_file = "mQTL_covar_10_methy_PC.txt",
    cis_outfile = "cis_mQTL_robinson.txt",
    snp_pos = "snp_pos.txt",
    probe_pos = "probe_pos.txt",
    model = "modelLINEAR"
  )
)
print(argv)
use_model <- get(argv$model)
# @TODO fill in functionality to have other input files and output strings etc
data_dir <- argv$data_dir
SNP_fname <- paste0(data_dir, argv$SNP_fname)
methylation_fname <- paste0(data_dir, argv$methylation_fname)
cis_outfile <- paste0(data_dir, argv$cis_outfile)
pv_out_threshold <- 1.0 # @TODO check this
error_cov <- numeric()
cis_dist <- 75000

methylation <- SlicedData$new()
methylation$fileDelimiter <- "\t"
methylation$fileOmitCharacters <- "NA"
methylation$fileSkipRows <- 1
methylation$fileSkipColumns <- 1
methylation$fileSliceSize <- 10000 # 2000 methylation at once
methylation$LoadFile(methylation_fname)
probe_pos <- read.delim(paste0(data_dir, argv$probe_pos), sep = "\t")

covariates <- SlicedData$new()
covariates$fileDelimiter <- " "
covariates$fileOmitCharacters <- "NA"
covariates$fileSkipRows <- 1
covariates$fileSkipColumns <- 1
covariates$LoadFile(paste0(data_dir, argv$cov_file))

snps <- SlicedData$new()
snps$fileDelimiter <- "\t"
snps$fileOmitCharacters <- "NA"
snps$fileSkipRows <- 1
snps$fileSkipColumns <- 1
snps$fileSliceSize <- 10000 # 10000 snps at once
snps$LoadFile(SNP_fname)
snp_pos <- read.delim(paste0(data_dir, argv$snp_pos), sep = "\t")
print(head(snp_pos))
print(head(probe_pos))
me <- Matrix_eQTL_main(
  snps = snps,
  gene = methylation,
  cvrt = covariates,
  output_file_name.cis = cis_outfile,
  pvOutputThreshold = 0,
  pvOutputThreshold.cis = pv_out_threshold,
  cisDist = cis_dist,
  snpspos = snp_pos,
  genepos = probe_pos,
  useModel = use_model,
  errorCovariance = error_cov,
  verbose = TRUE,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)
cat("Analysis done in: ", me$time.in.sec, " seconds", "\n")
cat("Detected local eQTLs:", "\n")
nrow(me$cis)
