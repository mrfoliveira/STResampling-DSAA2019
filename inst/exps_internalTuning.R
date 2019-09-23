# PARALLELIZATION
NCORES <- 2
NUM_SPLITS <- NCORES
NUM_THREADS <- 1
library(doParallel)
cat(paste("\nUsing", NCORES, "cores and up to", NCORES*NUM_THREADS, "ranger threads\n\n"))
registerDoParallel(cores=NCORES)

# LOADING LEARNING MODEL FUNCTIONS

library(ranger)
library(earth)
library(rpart)

# FILE PATHS

DATA_PATH <- "./extdata/"
UTILS_PATH <- "../R" # if package not installed
RESULTS_PATH <- "./"
if(!dir.exists(RESULTS_PATH)) dir.create(RESULTS_PATH)

# LOADING PACKAGE CODE

if(!("STResamplingDSAA") %in% installed.packages()){
  tosource <- list.files(UTILS_PATH, full.names = TRUE)
  for(f in tosource) source(f)
}else{
  library(STResamplingDSAA)
}

# LOADING DATA 

cat("\nLoading data sets...\n")

load(paste0(DATA_PATH, "dfs.Rdata"))
load(paste0(DATA_PATH, "inds_df.Rdata"))

# PARAMETRIZATION

models <- c("rpart", "earth", "ranger")
cpercs <- list()
cpercs[["stunder"]] <- seq(0.1, 0.9, 0.1)
cpercs[["stover"]] <- c(seq(0.1, 1, 0.1), 2)
alphas <- c(0, 0.25, 0.5, 0.75, 1)

resample.grid <- list()
resample.grid$stunder <- expand.grid(alpha=alphas, C.perc=cpercs$stunder)
resample.grid$stover <- expand.grid(alpha=alphas, C.perc=cpercs$stover)
resample.grid$under <- data.frame(C.perc=cpercs$stunder)
resample.grid$over <- data.frame(C.perc=cpercs$stover)

EST_PARS <- list(nfolds = 10, 
                 window = "growing", 
                 fold.alloc.proc = "Tblock_SPall", 
                 removeSP = FALSE, 
                 time="time", 
                 site_id="station",
                 .keepTrain = TRUE,
                 .parallel = FALSE)

BASE_EST_PARS <- EST_PARS
BASE_EST_PARS$.parallel <- TRUE

INT_EST_PARS <- list(nfolds = 10, 
                     fold.alloc.proc = "Tblock_SPall", 
                     time="time", 
                     site_id="station",
                     .keepTrain = TRUE,
                     .parallel = TRUE)

THR_REL <- 0.9

EVAL_PARS <- list(eval.function = eval_stats, 
                  cf=1.5, thr=THR_REL, beta=1,
                  .keptTrain = TRUE)

IN_EVAL_PARS <- list(eval.function = eval_stats, 
                  cf=1.5, thr=THR_REL, beta=1)

SEED <- 1234 #rep1

inds_df <- inds_df[c('MESApol', 'NCDCPprec', 'TCEQOozone', 'TCEQTtemp', 'TCEQWwind', 
                     'RURALpm10', 'BEIJno', 'BEIJpm10', 'BEIJwind', 'BEIJpm25')]

# RUN EXPERIMENTS

res <- list()
for(m in models){
  res[[m]] <- list()
  
  for(d in 1:length(inds_df)){
    
    WF_PARS <- list(model=m, min_train=2, handleNAs="centralImput", nORp = 0.2)
    if(m=="ranger")
      WF_PARS <- c(WF_PARS, list(num.trees = 250, num.threads=NUM_THREADS, verbose=FALSE))
    
    RS_PARS <- list(thr.rel=THR_REL)
    
    dfnm <- names(inds_df)[d]
    cat(paste("\n\nTesting data", dfnm, "with", m,"\n"))
    
    ind_df <- as.data.frame(inds_df[[dfnm]]$df)
    stations <- data_list[[dfnm]]$stations
    
    res[[m]][[dfnm]] <- list()
    
    cat("\nBaseline results...\n")

    
    res[[m]][[dfnm]][["baseline"]] <- estimates(ind_df, form=value~., 
                                            estimator="prequential_eval",
                                            est.pars = BASE_EST_PARS, 
                                            workflow = "simple_workflow", 
                                            wf.pars = WF_PARS, 
                                            evaluator = "evaluate", 
                                            eval.pars = EVAL_PARS, 
                                            seed=SEED)
    
    
    for(rsfun in c("under", "over", "stunder", "stover")){
      
      cat(paste("Internal validation for", rsfun,"\n"))
      
      if(rsfun == "under"){
        WF_PARS <- c(WF_PARS, 
                     list(internal.est = "kf_xval", 
                          internal.est.pars = INT_EST_PARS,
                          internal.evaluator = "int_util_evaluate", 
                          internal.eval.pars = IN_EVAL_PARS,
                          metrics = c("F1.u", "rmse_phi"), 
                          metrics.max = c(TRUE, FALSE),
                          stat="MEAN")) 
      }
      if(rsfun=="stunder"){
        # getting timestamps with right format
        if(grepl("BEIJ", dfnm)) ind_df$time <- lubridate::ymd_hms(ind_df$time)
        else ind_df$time <- as.Date(ind_df$time)
        
        RS_PARS <- c(RS_PARS, list(sites_sf=stations, type="add"))
      }
      
      WF_PARS$resample <- rsfun
      WF_PARS$resample.pars <- RS_PARS 
      WF_PARS$resample.grid <- resample.grid[[rsfun]]
      
      try( res[[m]][[dfnm]][[rsfun]] <- estimates(ind_df, form=value~., 
                                                  estimator="prequential_eval",
                                                  est.pars = EST_PARS, 
                                                  workflow = "internal_workflow", 
                                                  wf.pars = WF_PARS, 
                                                  evaluator = "evaluate", 
                                                  eval.pars = EVAL_PARS, 
                                                  seed=SEED)  
           )  
           
    }
    
    cat("\nSaving results...\n")
    save(res, file=paste0(RESULTS_PATH,"res_internalTuning.Rdata"))
  }
  
}
cat("\n\nDone!\n\n")


# PROCESS AND SAVE RESULTS

library(dplyr)

sumRes_int <- lapply(res, function(x) lapply(x, function(y) lapply(y, function(z) z$evalRes)))

sumResTab_int <- bind_rows(lapply(res, function(x) 
  bind_rows(lapply(x, function(y) 
    bind_rows(lapply(y, function(z) cbind(data.frame(fold=1:nrow(z$evalRes)), 
      as.data.frame(z$evalRes))), .id="sampling")), 
  .id="data")),
.id="model")

chosen <- bind_rows(
  lapply(res, function(m) dplyr::bind_rows(
    lapply(m, function(d) dplyr::bind_rows(
      lapply(d, function(y) dplyr::bind_rows(
        lapply(y$rawRes, function(x) x$resample_chosen))),
      .id="sampling")),
    .id="data")),
  .id="model") %>% 
select(model, data, sampling, C.perc, alpha) %>% 
rename(cperc=C.perc)

int_grid_res <- bind_rows(
  lapply(res, function(m) dplyr::bind_rows(
    lapply(m, function(d) dplyr::bind_rows(
      lapply(d, function(y) dplyr::bind_rows(
        lapply(y$rawRes, function(x) x$resample_grid))),
      .id="sampling")),
    .id="data")),
  .id="model") 

save(chosen, int_grid_res, sumRes_int, sumResTab_int, file=paste0(RESULTS_PATH, "sumRes_internalTuning.Rdata"))
