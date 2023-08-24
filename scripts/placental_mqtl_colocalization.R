library(tidyverse)
library(data.table)
library(coloc)
library(bigsnpr)
library(parallel)
library(glue)
argv <- commandArgs(trailingOnly = TRUE)
which_sex <- argv[[3]]
if (which_sex == "male_mqtl") {
  marginal_bonf <- fread("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/male_mcpg_bonf.txt.gz", key = "SNP")[p < 0.05]
} else if (which_sex == "female_mqtl") {
  marginal_bonf <- fread("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/female_mcpg_bonf.txt.gz", key = "SNP")[p < 0.05]
}else{
  marginal_bonf <- fread("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/marginal_mcpg_bonf.txt.gz", key = "SNP")[p < 0.05]
}
rds <- snp_readBed2("/arc/project/st-dennisjk-1/shared/data/1000G_EUR_ldsc_data/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL.bed", backingfile = tempfile())

reference <- snp_attach(rds)

compute_coloc <- function(SNP, mqtl, gwas, method = "coloc", type = "quant", s = NA) { # SNP must be in SNP column of mqtl and gwas
  D1 <- list(
    beta = mqtl$b,
    varbeta = mqtl$SE^2,
    snp = mqtl$SNP,
    position = mqtl$BP,
    N = 400,
    MAF = mqtl$Freq,
    type = "quant"
  )
  if (is.na(s)) {
    D2 <- list(
      pvalues = pnorm(-abs(gwas$Z)) * 2,
      z = gwas$Z,
      snp = gwas$SNP,
      MAF = snp_MAF(reference$genotypes, ind.col = match(gwas$SNP, reference$map$marker.ID)),
      N = min(gwas$N),
      type = type
    )
  } else {
    D2 <- list(
      pvalues = pnorm(-abs(gwas$Z)) * 2,
      z = gwas$Z,
      snp = gwas$SNP,
      MAF = snp_MAF(reference$genotypes, ind.col = match(gwas$SNP, reference$map$marker.ID)),
      N = min(gwas$N),
      type = type,
      s = s
    )
  }
  if (method == "susie") {
    LD <- snp_cor(reference$genotypes, ind.col = match(SNP, reference$map$marker.ID))^2
    LD <- as.matrix(LD)
    colnames(LD) <- SNP
    rownames(LD) <- SNP
    D1$LD <- LD
    D2$LD <- LD
    S1 <- runsusie(D1)
    S2 <- runsusie(D2)
    return(coloc.susie(S1, S2))
  } else if (method == "coloc") {
    return(coloc.abf(D1, D2))
  } else {
    return(NULL)
  }
}

if (argv[[1]] == "PGC") {
  sumstat_files <- c(
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/adhd_jul2017.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/anxiety.meta.full.fs.tbl.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/anxiety.meta.full.cc.tbl.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pgcAN2.2019-07.vcf.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/AUDIT_UKB_2018_AJP.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/Cannabis_ICC_23andmetop_UKB_het.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/iPSYCH-PGC_ASD_Nov2017.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/ocd_aug2017.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pgc_alcdep.eur_unrelated.aug2018_release.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pgc_bip_2018.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/PGC_UKB_depression_genome-wide.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/MDD2018_ex23andMe.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/PGC3_SCZ_wave3_public.v2.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/tag.cpd.tbl.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/TS_Oct2018.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pts_all_freeze2_overall.results.sumstats.gz"
  )
  trait_names <- c(
    "ADHD",
    "ANXFS",
    "ANX",
    "AN",
    "AUDIT",
    "CUD",
    "ASD",
    "OCD",
    "ALC",
    "BIP",
    "MDD",
    "MDD2018",
    "SCZ",
    "TAG_CPD",
    "TS",
    "PTSD"
  )
  sample_prev <- c(
    0.36,
    NA,
    0.33,
    0.23,
    NA,
    0.30,
    0.40,
    0.28,
    0.26,
    0.15,
    0.34,
    0.35,
    0.42,
    NA,
    0.34,
    0.15
  )
} else if (argv[[1]] == "neonatal") {
  sumstat_files <- c(
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/Pubertal_growth_PGF_PGM_combined.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG_Obesity_Meta_Analysis_1.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_hay_fever_irnt.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG-TotalGWG-Offspring.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/Pubertal_growth_PTF_PTM_combined.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_hay_fever_raw.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_asthma_raw.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_asthma_irnt.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_diabetes_raw.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_diabetes_irnt.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG_BMI_HapMap_DISCOVERY.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/Pubertal_growth_10F_12M_combined.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG_HC_DISCOVERY.v2.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/CHILD_ONSET_ASTHMA.20180501.allchr.assoc.GC.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/interpreggen.fetal.pe.meta.release.31jan2017.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/mat_all_chrALL_STERR_EU.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/ukbb_preeclampsia.gwas.imputed_v3.female.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/T1D.UCSC_META.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG-GWAS-BL.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/Fetal_BW_European_meta.NG2019.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/Fetal_Effect_European_meta_NG2019.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/Maternal_BW_European_meta.NG2019.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/hayfever_eczema_irnt.gwas.imputed_v3.both_sexes.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/hayfever_eczema_raw.gwas.imputed_v3.both_sexes.tsv.sumstats.gz"

  )
  trait_names <- c(
    "PGF_PGM",
    "EGG_OBESITY",
    "HAY_FEVER(irnt)",
    "GWG",
    "PTF_PTM",
    "HAY_FEVER(raw)",
    "AGE_ASTHMA(raw)",
    "AGE_ASTHMA(irnt)",
    "AGE_DIABETES(raw)",
    "AGE_DIABETES(irnt)",
    "EGG_BMI",
    "10F_12M",
    "EGG_GC",
    "CHILD_ASTHMA",
    "FETAL_PREECLAMPSIA",
    "MATERNAL_PREECLAMPSIA",
    "UKBB_PREECLAMPSIA",
    "Type 1 Diabetes",
    "EGG_BIRTH_LENGTH",
    "EGG_BIRTH_WEIGHT_FETAL",
    "EGG_BIRTH_WEIGHT_FETAL_EFFECT",
    "EGG_BIRTH_WEIGHT_MATERNAL",
    "HAY_FEVER_ECZEMA(irnt)",
    "HAY_FEVER_ECZEMA(raw)"
  )
  sample_prev <- c(
    NA,
    0.399,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    0.0288,
    0.00862,
    0.5,
    0.0108,
    0.0364,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA
  )
} else if (argv[[1]] == "male") {
  marginal_bonf <- fread("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/male_mcpg_bonf.txt.gz", key = "SNP")[p < 0.05]
  sumstat_files <- c(
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_asthma_irnt.gwas.imputed_v3.male.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_hay_fever_raw.gwas.imputed_v3.male.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_asthma_raw.gwas.imputed_v3.male.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG_TANNER_males.v2.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_diabetes_irnt.gwas.imputed_v3.male.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_diabetes_raw.gwas.imputed_v3.male.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/ukbb_preeclampsia.gwas.imputed_v3.male.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_hay_fever_irnt.gwas.imputed_v3.male.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/meta_STDERR_bip_eur_auto_M1_08_gcOFF_pgc.txt.gz.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/OCD_meta_male_auto_072416.gz.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/meta_STDERR_mdd_eur_auto_M1_08_gcOFF_pgc.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pgc_adhd_males.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/meta_STDERR_rmdd_eur_auto_M1_08_gcOFF_pgc.txt.gz.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pts_all_freeze2_males.results.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/meta_STDERR_scz_eur_auto_M1_08_gcOFF_pgc.txt.gz.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/ukbb_anxiety.gwas.imputed_v3.male.tsv.gz.fixed.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/NEW_META_PGC_iPSYCH_ASD_males.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/hayfever_eczema_irnt.gwas.imputed_v3.male.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/hayfever_eczema_raw.gwas.imputed_v3.male.tsv.sumstats.gz"
  )
  trait_names <- c(
    "AGE_ASTHMA_MALE(irnt)",
    "HAY_FEVER_MALE(raw)",
    "AGE_ASTHMA_MALE(raw)",
    "TANNER_STAGE_MALE",
    "AGE_DIABETES_MALE(irnt)",
    "AGE_DIABETES_MALE(raw)",
    "HAY_FEVER_MALE(irnt)",
    "BIP_MALE",
    "OCD_MALE",
    "MDD_MALE",
    "ADHD_MALE",
    "RMDD_MALE",
    "PTSD_MALE",
    "SCZ_MALE",
    "UKBB_ANXIETY_MALE",
    "ASD_MALE",
    "HAY_FEVER_ECZEMA_MALE(irnt)",
    "HAY_FEVER_ECZEMA_MALE(norm)"
  )
  sample_prev <- c(
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    0.351,
    0.309,
    0.309,
    0.441,
    0.245,
    0.127,
    0.517,
    0.0109,
    0.482,
    NA,
    NA
  )
} else if (argv[[1]] == "female") {
  marginal_bonf <- fread("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/female_mcpg_bonf.txt.gz", key = "SNP")[p < 0.05]
  sumstat_files <- c(
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_asthma_irnt.gwas.imputed_v3.female.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_hay_fever_raw.gwas.imputed_v3.female.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_asthma_raw.gwas.imputed_v3.female.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/EGG_TANNER_females.v2.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_diabetes_irnt.gwas.imputed_v3.female.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_diabetes_raw.gwas.imputed_v3.female.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/ukbb_preeclampsia.gwas.imputed_v3.female.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/age_hay_fever_irnt.gwas.imputed_v3.female.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/meta_STDERR_bip_eur_auto_F1_08_gcOFF_pgc.txt.gz.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/OCD_meta_female_auto_072416.gz.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/meta_STDERR_mdd_eur_auto_F1_08_gcOFF_pgc.txt.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pgc_adhd_females.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/meta_STDERR_rmdd_eur_auto_F1_08_gcOFF_pgc.txt.gz.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/pts_all_freeze2_females.results.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/meta_STDERR_scz_eur_auto_F1_08_gcOFF_pgc.txt.gz.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/ukbb_anxiety.gwas.imputed_v3.female.tsv.gz.fixed.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/pgc_formatted_sumstats/NEW_META_PGC_iPSYCH_ASD_females.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/mat_all_chrALL_STERR_EU.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/hayfever_eczema_irnt.gwas.imputed_v3.female.tsv.sumstats.gz",
    "/scratch/st-dennisjk-1/wcasazza/tmp_GWAS/neonatal_gwas/formatted/hayfever_eczema_raw.gwas.imputed_v3.female.tsv.sumstats.gz"
  )
  trait_names <- c(
    "AGE_ASTHMA_FEMALE(irnt)",
    "HAY_FEVER_FEMALE(raw)",
    "AGE_ASTHMA_FEMALE(raw)",
    "TANNER_STAGE_FEMALE",
    "AGE_DIABETES_FEMALE(irnt)",
    "AGE_DIABETES_FEMALE(raw)",
    "UKBB_PREECLAMPSIA_FEMALE",
    "HAY_FEVER_FEMALE(irnt)",
    "BIP_FEMALE",
    "OCD_FEMALE",
    "MDD_FEMALE",
    "ADHD_FEMALE",
    "RMDD_FEMALE",
    "PTSD_FEMALE",
    "SCZ_FEMALE",
    "UKBB_ANXIETY_FEMALE",
    "ASD_FEMALE",
    "MATERNAL_PREECLAMPSIA_FEMALE",
    "HAY_FEVER_ECZEMA_FEMALE(irnt)",
    "HAY_FEVER_ECZEMA_FEMALE(raw)"
  )
  sample_prev <- c(
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    0.0108,
    NA,
    0.43,
    0.261,
    0.477,
    0.233,
    0.402,
    0.283,
    0.37,
    0.0162,
    0.236,
    0.5,
    NA,
    NA
  )
}
result <- list()
i <- as.numeric(argv[[2]])
if (file.exists(glue("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/{trait_names[i]}_{which_sex}colocalization.txt"))) {
  quit()
}

gwas <- fread(sumstat_files[i])
tmp_marginal_bonf <- marginal_bonf[intersect(gwas$SNP, reference$map$marker.ID), on = "SNP", nomatch = 0]
tmp_marginal_bonf <- merge(tmp_marginal_bonf, gwas, by = "SNP")
print("Reading in elligble CpG sites")
eligible_cpg <- tmp_marginal_bonf[max(abs(Z)) > 5.45, .(Zmax = max(abs(Z)), minp = min(p)), by = "Probe"][Zmax > 5.45 & minp < 5e-8]$Probe
print(length(eligible_cpg))
if (length(eligible_cpg) > 0) {
  print("starting colocalization")
  test <- mclapply(
    eligible_cpg,
    function(probe) {
      mqtl <- tmp_marginal_bonf[probe, on = "Probe"]
      gwas_tmp <- gwas[mqtl$SNP, on = "SNP"]
      invisible(
        capture.output(res <- compute_coloc(
          mqtl$SNP,
          mqtl,
          gwas_tmp,
          method = "coloc",
          s = sample_prev[i],
          type = ifelse(is.na(sample_prev[i]), "quant", "cc")
        )$summary)
      )
      return(res)
    },
    mc.cores = 32
  )
  names(test) <- eligible_cpg
  result <- rbindlist(lapply(test, function(x) data.table(t(x))), idcol = "probe")

  fwrite(
    result,
    glue("/scratch/st-dennisjk-1/wcasazza/sex_specific_mQTL/data/{trait_names[i]}_{which_sex}colocalization.txt"),
    quote = F,
    sep = "\t",
    row.names = F
  )
}
