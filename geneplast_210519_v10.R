install.packages("BiocManager")
BiocManager::install("geneplast")
library(geneplast)

#load datasets
load("STRING10_2017.RData")


install.packages("readxl")
library(readxl)

cog_ids <-read_xlsx("OGs_list_for_R_vs_10.5.xlsx")


cogid_vec_apop <- cog_ids$apoptosis
cogid_vec_senesc <- cog_ids$senescence
cogid_vec_auto <- cog_ids$autophagy
cogid_vec_fanconi <- cog_ids$Fanconi_anemia_pathway
cogid_vec_nonhom <- cog_ids$`Non-homologous_end-joining`
cogid_vec_homrecomb <- cog_ids$Homologous_recombination
cogid_vec_mismatch <- cog_ids$Mismatch_repair
cogid_vec_baseexc <- cog_ids$Base_excision_repair
cogid_vec_dnarep <- cog_ids$DNA_replication
cogid_vec_nucleotideexc <- cog_ids$Nucleotide_excision_repair
cogid_vec_cellfate <- cog_ids$cell_fate
cogid_vec_genstab <- cog_ids$genome_maintenance
cogid_vec_total <- cog_ids$total

#create an object of class 'OGP'
ogp_apop <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_apop)
ogp_auto <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_auto)
ogp_baseexc <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_baseexc)
ogp_cellfate <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_cellfate)
ogp_dnarep <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_dnarep)
ogp_fanconi <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_fanconi)
ogp_homrecomb <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_homrecomb)
ogp_mismatch <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_mismatch)
ogp_nonhom <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_nonhom)
ogp_nucleotideexc <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_nucleotideexc)
ogp_genstab <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_genstab)
ogp_senesc <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_senesc)
ogp_total <- gplast.preprocess(cogdata=cogdata, sspids=sspids, cogids=cogid_vec_total)

## run the gplast function
ogp_apop <- gplast(ogp_apop)
ogp_auto <- gplast(ogp_auto)
ogp_baseexc <- gplast(ogp_baseexc)
ogp_cellfate <- gplast(ogp_cellfate)
ogp_dnarep <- gplast(ogp_dnarep)
ogp_fanconi <- gplast(ogp_fanconi)
ogp_homrecomb <- gplast(ogp_homrecomb)
ogp_mismatch <- gplast(ogp_mismatch)
ogp_nonhom <- gplast(ogp_nonhom)
ogp_nucleotideexc <- gplast(ogp_nucleotideexc)
ogp_genstab <- gplast(ogp_genstab)
ogp_senesc <- gplast(ogp_senesc)
ogp_total <- gplast(ogp_total)


res_apop <- gplast.get(ogp_apop,what="results")
res_auto <- gplast.get(ogp_auto,what="results")
res_baseexc <- gplast.get(ogp_baseexc,what="results")
res_cellfate <- gplast.get(ogp_cellfate,what="results")
res_dnarep <- gplast.get(ogp_dnarep,what="results")
res_fanconi <- gplast.get(ogp_fanconi,what="results")
res_homrecomb <- gplast.get(ogp_homrecomb,what="results")
res_mismatch <- gplast.get(ogp_mismatch,what="results")
res_nonhom <- gplast.get(ogp_nonhom,what="results")
res_nucleotideexc <- gplast.get(ogp_nucleotideexc,what="results")
res_genstab <- gplast.get(ogp_genstab,what="results")
res_senesc <- gplast.get(ogp_senesc,what="results")
res_total <- gplast.get(ogp_total,what="results")


write.csv(res_apop,"gplast_res_apop_10_diff.csv")
write.csv(res_auto, "gplast_res_auto_10_diff.csv")
write.csv(res_baseexc, "gplast_res_baseexc_10.csv")
write.csv(res_cellfate, "gplast_res_cellfate_10.csv")
write.csv(res_dnarep, "gplast_res_dnarep_10.csv")
write.csv(res_fanconi, "gplast_res_fanconi_10.csv")
write.csv(res_homrecomb, "gplast_res_homrecomb_10.csv")
write.csv(res_mismatch, "gplast_res_mismatch_10.csv")
write.csv(res_nonhom, "gplast_res_nonhom_10.csv")
write.csv(res_nucleotideexc, "gplast_res_nucleotideexc_11.csv")
write.csv(res_genstab, "gplast_res_genstab_10_diff.csv")
write.csv(res_senesc, "gplast_res_senesc_10_diff.csv")
write.csv(res_total, "gplast_res_total_10.csv")


#create an object of class 'OGR' 
ogr_apop <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_apop)
ogr_senesc <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_senesc)
ogr_auto <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_auto)
ogr_fanconi <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_fanconi)
ogr_nonhom <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_nonhom)
ogr_homrecomb <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_homrecomb)
ogr_mismatch <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_mismatch)
ogr_baseexc <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_baseexc)
ogr_dnarep <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_dnarep)
ogr_nucleotideexc <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_nucleotideexc)
ogr_cellfate <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_cellfate)
ogr_genstab <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_genstab)
ogr_total <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogid_vec_total)


ogr_apop <- groot(ogr_apop, nPermutations=10000)
ogr_senesc <- groot(ogr_senesc, nPermutations=10000)
ogr_auto <- groot(ogr_auto, nPermutations=10000)
ogr_fanconi <- groot(ogr_fanconi, nPermutations=10000)
ogr_nonhom <- groot(ogr_nonhom, nPermutations=10000)
ogr_homrecomb <- groot(ogr_homrecomb, nPermutations=10000)
ogr_mismatch <- groot(ogr_mismatch, nPermutations=10000)
ogr_baseexc <- groot(ogr_baseexc, nPermutations=10000)
ogr_dnarep <- groot(ogr_dnarep, nPermutations=10000)
ogr_nucleotideexc <- groot(ogr_nucleotideexc, nPermutations=10000)
ogr_cellfate <- groot(ogr_cellfate, nPermutations=10000)
ogr_genstab <- groot(ogr_genstab, nPermutations=10000)
ogr_total <- groot(ogr_total, nPermutations=10000)


res_apop <- groot.get(ogr_apop,what="results")
res_auto <- groot.get(ogr_auto,what="results")
res_baseexc <- groot.get(ogr_baseexc,what="results")
res_cellfate <- groot.get(ogr_cellfate,what="results")
res_dnarep <- groot.get(ogr_dnarep,what="results")
res_fanconi <- groot.get(ogr_fanconi,what="results")
res_homrecomb <- groot.get(ogr_homrecomb,what="results")
res_mismatch <- groot.get(ogr_mismatch,what="results")
res_nonhom <- groot.get(ogr_nonhom,what="results")
res_nucleotideexc <- groot.get(ogr_nucleotideexc,what="results")
res_genstab <- groot.get(ogr_genstab,what="results")
res_senesc <- groot.get(ogr_senesc,what="results")
res_total <- groot.get(ogr_total,what="results")


write.csv(res_apop,"groot_res_apop_10.csv")
write.csv(res_auto, "groot_res_auto_10.csv")
write.csv(res_baseexc, "groot_res_baseexc_10.csv")
write.csv(res_cellfate, "groot_res_cellfate_10.csv")
write.csv(res_dnarep, "groot_res_dnarep_10.csv")
write.csv(res_fanconi, "groot_res_fanconi_10.csv")
write.csv(res_homrecomb, "groot_res_homrecomb_10.csv")
write.csv(res_mismatch, "groot_res_mismatch_10.csv")
write.csv(res_nonhom, "groot_res_nonhom_10.csv")
write.csv(res_nucleotideexc, "groot_res_nucleotideexc_10.csv")
write.csv(res_genstab, "groot_res_genstab_10.csv")
write.csv(res_senesc, "groot_res_senesc_10.csv")
write.csv(res_total, "groot_res_total_10.csv")


lapply(cogid_vec_total, groot.plot, ogr=ogr_total, width = 10, height = 24, cex.lab = 0.52, cex.nodes = 1)





