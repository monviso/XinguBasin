#!/usr/bin/env Rscript

# STEP 01b — Fit the prescreen-selected hierarchy and assign all eligible shots
# NMF identifies relative RH-profile shape; within-family k-means separates
# absolute RH98 and negative depth. Annual inputs are already quality filtered.

options(stringsAsFactors = FALSE, warn = 1)

cfg <- list(
  input_files = file.path(
    "E:/Xingu_rev/data",
    sprintf("GEDI_RFLD_sampling_frame_%d.parquet", 2019:2023)
  ),
  prescreen_dir = "E:/Xingu_rev/data/step01a_classification_prescreen",
  output_dir = "E:/Xingu_rev/data/step01b_final_classification",
  seed = 42L,
  overwrite = FALSE,
  # Rebuild only Phase 3 while reusing the expensive sample and fitted model.
  overwrite_robustness = FALSE,
  
  # ---------------- EXACT STRATIFICATION DEFINITION ---------------- #
  training_target_n = 500000L,
  years = 2019:2023,
  mapbiomas_forest_classes = c(3L, 4L, 5L, 6L),
  d1_breaks_m = c(30, 100, 200, 500, 1000, 3000, Inf),
  d1_labels = c("30_100", "100_200", "200_500", "500_1000",
                "1000_3000", "gt_3000"),
  d2_breaks_m = c(0, 100, 200, 500, 1000, 3000, Inf),
  d2_labels = c("0_100", "100_200", "200_500", "500_1000",
                "1000_3000", "gt_3000"),
  patch_area_breaks_ha = c(0, 10, 100, 1000, 10000, Inf),
  patch_area_labels = c("lt_10", "10_100", "100_1000",
                        "1000_10000", "gt_10000"),
  compactness_breaks = c(0, .2, .4, .6, .8, 1, Inf),
  compactness_labels = c("0_.2", ".2_.4", ".4_.6", ".6_.8", ".8_1", "gt_1"),
  spatial_cell_size_m = 10000,
  minimum_per_landscape_stratum = 20L,
  maximum_per_landscape_stratum = 500L,
  allocation_exponent = 0.5,
  adjacent_lulc_groups = list(
    wetland = 11L, grassland = 12L, pasture = 15L,
    agriculture = c(18L, 19L, 20L, 39L, 40L, 41L, 46L, 47L, 48L),
    silviculture_oil_palm = c(9L, 35L), urban = 24L, mining = 30L,
    water = c(31L, 32L, 33L),
    other = c(13L, 21L, 23L, 25L, 27L, 29L, 49L, 50L, 62L, 68L)
  ),
  
  # NMF and height subdivision
  rank_selection_n = 100000L,
  ranks = 1:10,
  fixed_nmf_rank = NA_integer_,
  rank_restarts = 3L,
  final_restarts = 5L,
  nmf_maxit = 200L,
  nmf_tol = 1e-5,
  k_candidates = 1:6,
  fixed_k_by_component = NULL,
  expected_final_classes = NA_integer_,
  # Lloyd avoids Hartigan-Wong Quick-Transfer failures caused by the many
  # identical/near-identical two-variable height-feature observations.
  kmeans_algorithm = "Lloyd",
  kmeans_nstart = 50L,
  kmeans_iter_max = 500L,
  pairwise_alpha = .05,
  
  # Robustness (reduce robustness_n/reps for a quick test)
  robustness_reps = 3L,
  robustness_fraction = .80,
  robustness_n = 100000L,
  rcppml_threads = max(1L, parallel::detectCores(logical = FALSE) - 1L)
)

packages <- c("arrow", "data.table", "dplyr", "ggplot2", "Matrix", "RcppML",
              "inflection", "scales", "clue")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages({
  library(arrow); library(data.table); library(dplyr); library(ggplot2)
  library(Matrix); library(RcppML)
})
RcppML::setRcppMLthreads(cfg$rcppml_threads)
set.seed(cfg$seed)

selection_file<-file.path(cfg$prescreen_dir,"selection","selected_hierarchy.csv")
if(!file.exists(selection_file)) stop("Run Step 01a first; missing: ",selection_file)
selection<-fread(selection_file)
required_selection<-c("selected_nmf_rank","total_classes")
if(nrow(selection)!=1L||length(setdiff(required_selection,names(selection))))
  stop("Invalid Step 01a selected_hierarchy.csv")
cfg$fixed_nmf_rank<-as.integer(selection$selected_nmf_rank)
k_columns<-grep("^k_component_[0-9]+$",names(selection),value=TRUE)
k_components<-as.integer(sub("k_component_","",k_columns))
cfg$fixed_k_by_component<-setNames(
  as.integer(unlist(selection[1,..k_columns],use.names=FALSE)),
  as.character(k_components)
)
cfg$expected_final_classes<-as.integer(selection$total_classes)
if(!length(k_columns)||sum(cfg$fixed_k_by_component)!=cfg$expected_final_classes)
  stop("Step 01a hierarchy is internally inconsistent")

dirs <- list(
  checkpoints=file.path(cfg$output_dir,"checkpoints"),
  audit=file.path(cfg$output_dir,"audit"), models=file.path(cfg$output_dir,"models"),
  figures=file.path(cfg$output_dir,"figures"),
  robustness=file.path(cfg$output_dir,"robustness"),
  assignments=file.path(cfg$output_dir,"assignments_dataset")
)
invisible(lapply(c(cfg$output_dir, unlist(dirs)), dir.create,
                 recursive=TRUE, showWarnings=FALSE))
stopifnot(all(file.exists(cfg$input_files)))

rh_names <- paste0("rh",0:98)
required <- c("sample_id","shot_number","year","latitude","longitude",
              "mapbiomas_class","external_edge_distance_m",
              "nearest_other_forest_distance_m","nearest_opposite_lulc_class",
              "focal_patch_area_ha","focal_patch_compactness","tmf_class",rh_names)
ds <- arrow::open_dataset(cfg$input_files, format="parquet")
absent <- setdiff(required, names(ds))
if (length(absent)) stop("Missing input columns: ",paste(absent,collapse=", "))

logmsg <- function(...) {cat(sprintf("[%s] ",format(Sys.time(),"%F %T")),...,"\n",sep="");flush.console()}
cut_safe <- function(x,b,l) cut(x,b,labels=l,right=FALSE,include.lowest=TRUE)
adjacent_group <- function(x) {
  out <- rep("unclassified",length(x)); x <- as.integer(x)
  for (nm in names(cfg$adjacent_lulc_groups)) out[x %in% cfg$adjacent_lulc_groups[[nm]]] <- nm
  out[is.na(x)] <- "none_or_unavailable"; out
}
shape_normalize <- function(rh) {
  lo <- apply(rh,1,min,na.rm=TRUE); hi <- rh[,99L]; span <- hi-lo
  if (any(!is.finite(span)|span<=0)) stop("Invalid RH range in retained profiles")
  z <- (rh-lo)/span
  z[z<0 & z> -1e-6] <- 0; z[z>1 & z<1+1e-6] <- 1
  z <- pmin(pmax(z,0),1); colnames(z) <- paste0("rhs",0:98); z
}
as_sparse <- function(x) methods::as(Matrix::Matrix(x,sparse=TRUE),"dgCMatrix")
rowmax <- function(x) max.col(x,ties.method="first")
ari <- function(x,y) {
  z<-table(x,y); c2<-function(a)a*(a-1)/2; n<-sum(z)
  a<-sum(c2(z)); b<-sum(c2(rowSums(z))); c<-sum(c2(colSums(z)))
  e<-b*c/c2(n); m<-(b+c)/2; if(m==e) 1 else (a-e)/(m-e)
}
nmi <- function(x,y) {
  z<-table(x,y)/length(x); px<-rowSums(z); py<-colSums(z); ij<-which(z>0,arr.ind=TRUE)
  mi<-sum(z[ij]*log(z[ij]/(px[ij[,1]]*py[ij[,2]]))); hx<--sum(px[px>0]*log(px[px>0])); hy<--sum(py[py>0]*log(py[py>0]))
  if(hx==0&&hy==0) 1 else if(hx==0||hy==0) 0 else mi/sqrt(hx*hy)
}
choose_knee <- function(x,y) as.integer(x[which.min(abs(x-inflection::uik(x,y)))])

# Match a refitted NMF's arbitrary component labels to the reference model.
# This is essential before applying the fixed 3,2,2 within-family hierarchy.
match_nmf_components <- function(reference_w, candidate_w) {
  if (ncol(reference_w) != ncol(candidate_w)) stop("NMF ranks differ")
  rank<-ncol(reference_w)
  correlation <- cor(reference_w, candidate_w)
  correlation[!is.finite(correlation)]<- -1
  reference_to_candidate<-as.integer(clue::solve_LSAP(correlation-min(correlation),maximum=TRUE))
  candidate_to_reference <- integer(rank)
  candidate_to_reference[reference_to_candidate] <- seq_len(rank)
  list(
    reference_to_candidate=reference_to_candidate,
    candidate_to_reference=candidate_to_reference,
    mean_matched_correlation=mean(
      correlation[cbind(seq_len(rank),reference_to_candidate)],na.rm=TRUE
    )
  )
}

allocate_quota <- function(pop) {
  pop[, quota := pmin(N, cfg$minimum_per_landscape_stratum)]
  target <- min(cfg$training_target_n,sum(pop$N)); left <- target-sum(pop$quota)
  while(left>0) {
    pop[,capacity:=pmax(0,pmin(N,cfg$maximum_per_landscape_stratum)-quota)]
    active<-which(pop$capacity>0); if(!length(active)) break
    w<-pop$N[active]^cfg$allocation_exponent
    add<-pmin(pop$capacity[active],pmax(1L,floor(left*w/sum(w))))
    if(sum(add)>left) add<-add*0L+as.integer(seq_along(add)<=left)
    pop$quota[active]<-pop$quota[active]+add; left<-target-sum(pop$quota)
  }
  pop[,capacity:=NULL]; pop
}

# PHASE 1: build the landscape frame and spatially dispersed discovery sample.
sample_file <- file.path(cfg$prescreen_dir,"checkpoints","01_stratified_sample.rds")
if(!file.exists(sample_file)) stop("Missing Step 01a stratified-sample checkpoint")
if (FALSE) {
  logmsg("Collecting landscape variables for stratification")
  meta <- ds %>% select(all_of(setdiff(required,rh_names))) %>% collect() %>% as.data.table()
  meta <- meta[year %in% cfg$years & mapbiomas_class %in% cfg$mapbiomas_forest_classes]
  meta[,d1_bin:=cut_safe(external_edge_distance_m,cfg$d1_breaks_m,cfg$d1_labels)]
  meta[,d2_bin:=cut_safe(nearest_other_forest_distance_m,cfg$d2_breaks_m,cfg$d2_labels)]
  meta[,patch_area_bin:=cut_safe(focal_patch_area_ha,cfg$patch_area_breaks_ha,cfg$patch_area_labels)]
  meta[,compactness_bin:=cut_safe(focal_patch_compactness,cfg$compactness_breaks,cfg$compactness_labels)]
  meta[,adjacent_group:=adjacent_group(nearest_opposite_lulc_class)]
  mean_lat <- mean(meta$latitude,na.rm=TRUE)
  meta[,spatial_cell:=paste(floor(longitude*111320*cos(mean_lat*pi/180)/cfg$spatial_cell_size_m),
                            floor(latitude*110540/cfg$spatial_cell_size_m),sep="_")]
  meta[,landscape_stratum:=paste(year,mapbiomas_class,d1_bin,adjacent_group,
                                 patch_area_bin,sep="__")]
  pop<-meta[,.(N=.N),by=landscape_stratum]; pop<-allocate_quota(pop)
  set.seed(cfg$seed); meta[,random_order:=runif(.N)]
  first<-meta[order(random_order),.SD[!duplicated(spatial_cell)],by=landscape_stratum]
  first<-merge(first,pop[,.(landscape_stratum,quota)],by="landscape_stratum")
  first<-first[,head(.SD,unique(quota)),by=landscape_stratum]
  chosen<-first$sample_id
  first_n<-first[,.(first_n=.N),by=landscape_stratum]
  need<-merge(pop[,.(landscape_stratum,quota)],first_n,by="landscape_stratum",all.x=TRUE)
  need[is.na(first_n),first_n:=0L];need[,need:=quota-first_n]
  candidates<-merge(meta[!sample_id %in% chosen],need[,.(landscape_stratum,need)],by="landscape_stratum")
  fill<-candidates[need>0,head(.SD,unique(need)),by=landscape_stratum]
  selected<-unique(c(chosen,fill$sample_id))
  design<-meta[sample_id %in% selected]
  sampled_n<-design[,.(sampled_n=.N),by=landscape_stratum]
  design<-merge(design,pop[,.(landscape_stratum,population_n=N)],by="landscape_stratum")
  design<-merge(design,sampled_n,by="landscape_stratum")
  design[,selection_probability:=sampled_n/population_n]
  design[,sampling_weight:=population_n/sampled_n]
  fwrite(pop,file.path(dirs$audit,"stratification_population_and_quota.csv"))
  fwrite(design,file.path(dirs$audit,"stratification_selected_points.csv"))
  marginal_variables<-c("year","mapbiomas_class","d1_bin","d2_bin",
                        "adjacent_group","patch_area_bin","compactness_bin",
                        "tmf_class")
  marginal<-rbindlist(lapply(marginal_variables,function(v){
    population<-meta[,.(population_n=.N),by=v];setnames(population,v,"level")
    sampled<-design[,.(sampled_n=.N),by=v];setnames(sampled,v,"level")
    ans<-merge(population,sampled,by="level",all=TRUE)
    ans[is.na(sampled_n),sampled_n:=0L]
    ans[,`:=`(variable=v,population_proportion=population_n/sum(population_n),
              sample_proportion=sampled_n/sum(sampled_n))]
    ans[]
  }),use.names=TRUE,fill=TRUE)
  fwrite(marginal,file.path(dirs$audit,"stratification_marginal_balance.csv"))
  spatial_audit<-rbind(
    meta[,.(dataset="population",n=.N,unique_cells=uniqueN(spatial_cell)),by=year],
    design[,.(dataset="training_sample",n=.N,unique_cells=uniqueN(spatial_cell)),by=year]
  )
  fwrite(spatial_audit,file.path(dirs$audit,"stratification_spatial_coverage.csv"))
  capture.output(str(cfg),file=file.path(dirs$audit,"stratification_configuration.txt"))
  saveRDS(design,sample_file,compress=FALSE)
  rm(meta,first,fill,pop);gc()
}

# Read RH only for selected IDs.
training_file <- file.path(cfg$prescreen_dir,"checkpoints","02_training_profiles.rds")
if(!file.exists(training_file)) stop("Missing Step 01a RH-profile checkpoint")
if (FALSE) {
  design<-readRDS(sample_file); ids<-design$sample_id
  train<-ds %>% filter(sample_id %in% ids) %>% select(sample_id,all_of(rh_names)) %>% collect() %>% as.data.table()
  setkey(train,sample_id); setkey(design,sample_id); train<-design[train]
  rh<-as.matrix(train[,..rh_names]); storage.mode(rh)<-"double"
  rhs<-shape_normalize(rh)
  features<-scale(cbind(rh98=rh[,99],negative_depth=pmax(0,-rh[,1])))
  feature_center<-attr(features,"scaled:center")
  feature_scale<-attr(features,"scaled:scale");feature_scale[!is.finite(feature_scale)|feature_scale==0]<-1
  features<-sweep(cbind(rh98=rh[,99],negative_depth=pmax(0,-rh[,1])),2,feature_center,"-")
  features<-sweep(features,2,feature_scale,"/")
  saveRDS(list(meta=train[,setdiff(names(train),rh_names),with=FALSE],raw_rh=rh,rhs=rhs,
               height_features=features,feature_center=feature_center,
               feature_scale=feature_scale),training_file,compress=FALSE)
  rm(train,rh,rhs,features,design);gc()
}

fit_best_nmf <- function(A,k,restarts,offset=0L) {
  best<-NULL; bestm<-Inf; runs<-list()
  for(j in seq_len(restarts)) {
    f<-RcppML::nmf(A,k=k,tol=cfg$nmf_tol,maxit=cfg$nmf_maxit,seed=cfg$seed+offset+j,verbose=FALSE)
    m<-RcppML::mse(A,f$w,f$d,f$h); runs[[j]]<-data.table(restart=j,mse=m,iterations=f$iter)
    if(m<bestm){best<-f;bestm<-m}
  }; list(model=best,mse=bestm,runs=rbindlist(runs))
}
fit_kmeans <- function(component,x,offset=0L) {
  fits<-list(); rows<-list(); cluster<-integer(nrow(x))
  for(comp in sort(unique(component))) {
    ii<-which(component==comp); candidates<-cfg$k_candidates[cfg$k_candidates<=nrow(unique(as.data.table(x[ii,,drop=FALSE])))]
    local<-list(); w<-numeric(length(candidates))
    for(j in seq_along(candidates)){k<-candidates[j];set.seed(cfg$seed+offset+comp*100+k);local[[as.character(k)]]<-kmeans(x[ii,,drop=FALSE],k,nstart=cfg$kmeans_nstart,iter.max=cfg$kmeans_iter_max,algorithm=cfg$kmeans_algorithm);w[j]<-local[[as.character(k)]]$tot.withinss}
    component_name<-as.character(comp)
    if(!component_name %in% names(cfg$fixed_k_by_component)) stop("No fixed k for NMF component ",comp)
    k<-as.integer(cfg$fixed_k_by_component[[component_name]])
    if(!k %in% candidates) stop("Fixed k=",k," unavailable for NMF component ",comp)
    f<-local[[as.character(k)]]; ord<-order(f$centers[,1],f$centers[,2]); rel<-integer(k);rel[ord]<-seq_len(k);cluster[ii]<-rel[f$cluster]
    fits[[as.character(comp)]]<-list(k=k,centers=f$centers[ord,,drop=FALSE])
    rows[[as.character(comp)]]<-data.table(nmf_component=comp,k=candidates,wcss=w,selected=candidates==k)
  };list(cluster=cluster,fits=fits,wcss=rbindlist(rows))
}

# PHASE 2: rank, final NMF and within-shape height clustering.
model_file<-file.path(dirs$models,"nmf_shape_height_kmeans_model.rds")
checkpoint_model<-file.path(dirs$checkpoints,"03_model_and_training_classes.rds")
if (!file.exists(checkpoint_model) || cfg$overwrite) {
  o<-readRDS(training_file)
  rank_runs_path<-file.path(cfg$prescreen_dir,"models","nmf_rank_restarts.csv")
  rank_selection_path<-file.path(cfg$prescreen_dir,"models","nmf_rank_selection.csv")
  if(!file.exists(rank_runs_path)||!file.exists(rank_selection_path))
    stop("Step 01a rank diagnostics are incomplete")
  rr<-fread(rank_runs_path);rs<-fread(rank_selection_path)
  selected_rank<-as.integer(cfg$fixed_nmf_rank)
  if(!selected_rank %in% rs$rank) stop("Selected NMF rank is absent from Step 01a diagnostics")
  rs[,selected:=rank==selected_rank]
  fwrite(rr,file.path(dirs$models,"nmf_rank_restarts.csv"));fwrite(rs,file.path(dirs$models,"nmf_rank_selection.csv"))
  final<-fit_best_nmf(as_sparse(t(o$rhs)),selected_rank,cfg$final_restarts,50000)
  prescreen_model_path<-file.path(cfg$prescreen_dir,"models","nmf_shape_height_kmeans_model.rds")
  if(!file.exists(prescreen_model_path)) stop("Missing Step 01a reference NMF model")
  prescreen_model<-readRDS(prescreen_model_path)
  component_match<-match_nmf_components(prescreen_model$nmf$w,final$model$w)
  p<-component_match$reference_to_candidate
  final$model$w<-final$model$w[,p,drop=FALSE]
  final$model$d<-final$model$d[p]
  final$model$h<-final$model$h[p,,drop=FALSE]
  comp<-rowmax(t(final$model$h));km<-fit_kmeans(comp,o$height_features,60000)
  map<-data.table(comp=comp,local=km$cluster,rh98=o$raw_rh[,99],neg=pmax(0,-o$raw_rh[,1]))[,.(n=.N,mean_rh98=mean(rh98),median_rh98=median(rh98),median_negative_depth=median(neg)),by=.(comp,local)][order(mean_rh98,median_rh98,median_negative_depth)]
  map[,architecture_class:=.I];map[,key:=paste(comp,local,sep="__")]
  if(nrow(map)!=cfg$expected_final_classes) stop("Fixed hierarchy produced ",nrow(map)," rather than ",cfg$expected_final_classes," classes")
  cls<-map$architecture_class[match(paste(comp,km$cluster,sep="__"),map$key)]
  bundle<-list(nmf=final$model,selected_rank=selected_rank,kmeans=km$fits,class_map=map,
               feature_center=o$feature_center,feature_scale=o$feature_scale,
               normalization="Per footprint: (RH-RHmin)/(RH98-RHmin)",
               hierarchy_definition=paste0(
                 "Step 01a-selected NMF rank ",selected_rank,"; component k: ",
                 paste(names(cfg$fixed_k_by_component),cfg$fixed_k_by_component,
                       sep="=",collapse=", ")
               ),
               prescreen_to_final_mean_basis_correlation=component_match$mean_matched_correlation,
               config=cfg)
  saveRDS(bundle,model_file);fwrite(km$wcss,file.path(dirs$models,"kmeans_wcss.csv"))
  saveRDS(c(o,list(component=comp,local_cluster=km$cluster,architecture_class=cls,model=bundle)),checkpoint_model,compress=FALSE)
  training_output<-cbind(o$meta,data.table(nmf_component=comp,
                                           within_component_cluster=km$cluster,architecture_class=cls,
                                           rh0_original=o$raw_rh[,1],rh98_original=o$raw_rh[,99],
                                           negative_depth=pmax(0,-o$raw_rh[,1])))
  arrow::write_parquet(training_output,file.path(cfg$output_dir,"nmf_training_sample.parquet"),compression="zstd")
  fwrite(training_output[,.(n=.N,median_rh98=median(rh98_original),
                            median_negative_depth=median(negative_depth)),by=.(architecture_class,nmf_component)],
         file.path(dirs$models,"training_class_summary.csv"))
}

# Figures: knees and normalized/original profile summaries.
o<-readRDS(checkpoint_model);rs<-fread(file.path(dirs$models,"nmf_rank_selection.csv"));kw<-fread(file.path(dirs$models,"kmeans_wcss.csv"))
ggsave(file.path(dirs$figures,"NMF_rank_knee.png"),ggplot(rs,aes(rank,mse))+geom_line()+geom_point()+geom_point(data=rs[selected==TRUE],colour="red",size=3)+geom_vline(data=rs[selected==TRUE],aes(xintercept=rank),linetype=2,colour="red")+theme_bw()+scale_x_continuous(breaks=cfg$ranks),width=7,height=5,dpi=300)
ggsave(file.path(dirs$figures,"kmeans_knees_by_NMF_family.png"),ggplot(kw,aes(k,wcss))+geom_line()+geom_point()+geom_point(data=kw[selected==TRUE],colour="red",size=3)+facet_wrap(~nmf_component,scales="free_y")+theme_bw(),width=10,height=5,dpi=300)

# Original-scale RH summaries requested for ecological interpretation.
rh_report<-data.table(
  architecture_class=o$architecture_class,
  rh25=o$raw_rh[,26L],rh50=o$raw_rh[,51L],
  rh75=o$raw_rh[,76L],rh98=o$raw_rh[,99L]
)[,.(
  n=.N,
  rh25_mean=mean(rh25),rh25_sd=sd(rh25),
  rh50_mean=mean(rh50),rh50_sd=sd(rh50),
  rh75_mean=mean(rh75),rh75_sd=sd(rh75),
  rh98_mean=mean(rh98),rh98_sd=sd(rh98)
),by=architecture_class][order(architecture_class)]
fwrite(rh_report,file.path(dirs$models,"final_class_RH25_RH50_RH75_RH98_mean_sd.csv"))

# Definitive class comparisons on original-scale RH quantiles. Class labels are
# already globally ordered by mean RH98: class 1 is shortest and class N tallest.
rh_test_data<-data.table(
  architecture_class=factor(o$architecture_class),
  rh25=o$raw_rh[,26L],rh50=o$raw_rh[,51L],
  rh75=o$raw_rh[,76L],rh98=o$raw_rh[,99L]
)
rh_metrics<-c("rh25","rh50","rh75","rh98")
rh_omnibus<-rbindlist(lapply(rh_metrics,function(metric){
  z<-rh_test_data[is.finite(get(metric))]
  test<-kruskal.test(z[[metric]],z$architecture_class)
  data.table(metric=metric,statistic=as.numeric(test$statistic),
             df=as.numeric(test$parameter),p_value=test$p.value,
             significant_0_05=test$p.value<cfg$pairwise_alpha)
}))
rh_pairs<-t(combn(levels(rh_test_data$architecture_class),2L))
rh_pairwise<-rbindlist(lapply(rh_metrics,function(metric){
  ans<-rbindlist(lapply(seq_len(nrow(rh_pairs)),function(j){
    a<-rh_pairs[j,1L];b<-rh_pairs[j,2L]
    x<-rh_test_data[architecture_class==a & is.finite(get(metric)),get(metric)]
    y<-rh_test_data[architecture_class==b & is.finite(get(metric)),get(metric)]
    test<-wilcox.test(x,y,exact=FALSE,conf.int=FALSE)
    pair_count<-as.double(length(x))*as.double(length(y))
    superiority<-as.numeric(test$statistic)/pair_count
    data.table(metric=metric,class_a=as.integer(a),class_b=as.integer(b),
               n_a=length(x),n_b=length(y),mean_a=mean(x),mean_b=mean(y),
               mean_difference=mean(x)-mean(y),median_a=median(x),median_b=median(y),
               p_value=test$p.value,rank_biserial=2*superiority-1)
  }))
  ans[,p_adjusted_BH:=p.adjust(p_value,method="BH")]
  ans[,significant_0_05:=p_adjusted_BH<cfg$pairwise_alpha]
  ans[]
}))
fwrite(rh_omnibus,file.path(dirs$models,"final_RH_quantile_Kruskal_Wallis.csv"))
fwrite(rh_pairwise,file.path(dirs$models,"final_RH_quantile_pairwise_Wilcoxon_BH.csv"))
profile_summary<-function(mat,g,scale_name){rbindlist(lapply(sort(unique(g)),function(gr){x<-mat[g==gr,,drop=FALSE];data.table(group=gr,percentile=0:98,mean=colMeans(x),median=apply(x,2,median),lo=apply(x,2,quantile,.1),hi=apply(x,2,quantile,.9),scale=scale_name)}))}
p1<-profile_summary(o$rhs,o$component,"Relative profile shape");p2<-profile_summary(o$raw_rh,o$architecture_class,"Original RH")
plotprof<-function(d,y,fn){ggsave(file.path(dirs$figures,fn),ggplot(d,aes(percentile,median,colour=factor(group),fill=factor(group)))+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.12,colour=NA)+geom_line(linewidth=.8)+geom_line(aes(y=mean),linetype=2,linewidth=.45)+theme_bw()+labs(x="RH percentile",y=y,colour="Group",fill="Group",caption="Solid: median; dashed: mean; ribbon: 10th–90th percentiles"),width=9,height=6,dpi=300)}
plotprof(p1,"Normalized relative height","NMF_shape_family_profiles.png");plotprof(p2,"Height relative to ground (m)","final_class_original_RH_profiles.png")

# PHASE 3: subsampling robustness on a bounded pool.
rob_file<-file.path(dirs$robustness,"subsample_stability.csv")
if(!file.exists(rob_file)||cfg$overwrite||cfg$overwrite_robustness){
  set.seed(cfg$seed+99999L)
  pool<-sample.int(nrow(o$rhs),min(cfg$robustness_n,nrow(o$rhs)))
  rows<-list()
  for(r in seq_len(cfg$robustness_reps)){
    set.seed(cfg$seed+100000+r)
    ii<-sample(pool,floor(length(pool)*cfg$robustness_fraction))
    f<-fit_best_nmf(as_sparse(t(o$rhs[ii,,drop=FALSE])),o$model$selected_rank,
                    cfg$final_restarts,100000+r)
    raw_unmatched<-rowmax(t(f$model$h))
    component_match<-match_nmf_components(o$model$nmf$w,f$model$w)
    raw<-component_match$candidate_to_reference[raw_unmatched]
    kmr<-fit_kmeans(raw,o$height_features[ii,,drop=FALSE],200000+r)
    rep_class<-paste(raw,kmr$cluster,sep="__")
    reference_class_key<-paste(o$component[ii],o$local_cluster[ii],sep="__")
    rows[[r]]<-data.table(
      replicate=r,n=length(ii),mse=f$mse,
      selected_final_class_count=uniqueN(rep_class),
      mean_matched_basis_correlation=component_match$mean_matched_correlation,
      reference_component_ARI=ari(o$component[ii],raw),
      reference_component_NMI=nmi(o$component[ii],raw),
      final_class_ARI=ari(reference_class_key,rep_class),
      final_class_NMI=nmi(reference_class_key,rep_class)
    )
  }
  fwrite(rbindlist(rows),rob_file)
}

# PHASE 4: project every eligible footprint and save restartable batch files.
assignment_success<-file.path(dirs$assignments,"_SUCCESS")
if(!file.exists(assignment_success)||cfg$overwrite){
  bundle<-readRDS(model_file)
  w<-sweep(bundle$nmf$w,2,bundle$nmf$d,"*")
  gram<-crossprod(w)
  assignment_cols<-unique(c("sample_id","shot_number","year","latitude","longitude",
                            "mapbiomas_class","tmf_class","external_edge_distance_m",
                            "nearest_other_forest_distance_m","nearest_opposite_lulc_class",
                            "focal_patch_area_ha","focal_patch_compactness",
                            "nearest_nonforest_patch_area_ha",rh_names))
  assignment_cols<-intersect(assignment_cols,names(ds))
  reader<-arrow::as_record_batch_reader(ds %>% select(all_of(assignment_cols)))
  batch_id<-0L;total<-0L
  repeat{
    batch<-reader$read_next_batch();if(is.null(batch))break
    batch_id<-batch_id+1L;d<-as.data.table(as.data.frame(batch))
    raw<-as.matrix(d[,..rh_names]);storage.mode(raw)<-"double";rhs<-shape_normalize(raw)
    h<-RcppML::nnls(gram,crossprod(w,t(rhs)));if(is.null(dim(h)))h<-matrix(h,nrow=ncol(w))
    component<-rowmax(t(h))
    xf<-cbind(rh98=raw[,99],negative_depth=pmax(0,-raw[,1]))
    xf<-sweep(xf,2,bundle$feature_center,"-");xf<-sweep(xf,2,bundle$feature_scale,"/")
    local<-integer(nrow(d))
    for(comp in sort(unique(component))){ii<-which(component==comp);cent<-bundle$kmeans[[as.character(comp)]]$centers;distance<-sapply(seq_len(nrow(cent)),function(j)rowSums((xf[ii,,drop=FALSE]-matrix(cent[j,],nrow=length(ii),ncol=2,byrow=TRUE))^2));if(is.null(dim(distance)))distance<-matrix(distance,ncol=1);local[ii]<-rowmax(-distance)}
    key<-paste(component,local,sep="__");cls<-bundle$class_map$architecture_class[match(key,bundle$class_map$key)]
    out<-d[,setdiff(names(d),rh_names),with=FALSE]
    out[,`:=`(rh0_original=raw[,1],rh98_original=raw[,99],negative_depth=pmax(0,-raw[,1]),nmf_component=component,within_component_cluster=local,architecture_class=as.integer(cls))]
    hp<-as.data.table(t(h));setnames(hp,paste0("nmf_score_",seq_len(ncol(hp))));out<-cbind(out,hp)
    arrow::write_parquet(out,file.path(dirs$assignments,sprintf("assignment_batch_%06d.parquet",batch_id)),compression="zstd")
    total<-total+nrow(out);if(batch_id==1L||batch_id%%20L==0L)logmsg("Assigned ",format(total,big.mark=",")," footprints")
    rm(d,raw,rhs,h,xf,out,hp);gc(FALSE)
  }
  writeLines(c(paste0("rows=",total),paste0("completed=",Sys.time())),assignment_success)
}


#Determining NMF-final classes trajectory
model <- readRDS(
  "E:/Xingu_rev/data/step01b_final_classification/models/nmf_shape_height_kmeans_model.rds"
)

class_mapping <- data.table::as.data.table(model$class_map)

class_mapping[
  order(architecture_class),
  .(
    architecture_class,
    nmf_component = comp,
    within_NMF_cluster = local,
    n,
    mean_rh98,
    median_rh98,
    median_negative_depth
  )
]

data.table::fwrite(
  class_mapping[order(architecture_class)],
  "E:/Xingu_rev/data/step01b_final_classification/models/final_class_NMF_mapping.csv"
)
writeLines(capture.output(sessionInfo()),file.path(cfg$output_dir,"sessionInfo.txt"))
logmsg("STEP 01b complete: definitive model, diagnostics and all-shot assignments written to ",cfg$output_dir)
