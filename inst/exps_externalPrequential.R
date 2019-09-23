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
cpercs[["stunder"]] <- c(0.2, 0.4, 0.6, 0.8, 0.95)
cpercs[["stover"]] <- c(0.5,1,2,3,4)

alphas <- c(0, 0.25, 0.5, 0.75, 1)

inds_df <- inds_df[c('MESApol', 'NCDCPprec', 'TCEQOozone', 'TCEQTtemp', 'TCEQWwind', 
                     'RURALpm10', 'BEIJno', 'BEIJpm10', 'BEIJwind', 'BEIJpm25')]


EST_PARS <- list(nfolds = 10, 
                 window = "growing", 
                 fold.alloc.proc = "Tblock_SPall", 
                 removeSP = FALSE, 
                 time="time", 
                 site_id="station",
                 .keepTrain = TRUE,
                 .parallel = TRUE)

THR_REL <- 0.9

EVAL_PARS <- list(eval.function = eval_stats, 
  cf=1.5, thr=THR_REL, beta=1,
  .keptTrain = TRUE)

SEED <- 1234 # rep1

# RUN EXPERIMENTS

res <- list()
for(m in models){
  res[[m]] <- list()

  WF_PARS <- list(model=m, min_train=2, handleNAs="centralImput", nORp = 0.2)

  if(m=="ranger")
    WF_PARS <- c(WF_PARS, list(num.trees = 250, num.threads=NUM_THREADS, verbose=FALSE))
  
  for(d in 1:length(inds_df)){
    
    dfnm <- names(inds_df)[d]
    cat(paste("\n\nTesting data", dfnm, "and model", m, "\n"))
    
    ind_df <- as.data.frame(inds_df[[dfnm]]$df)
    stations <- data_list[[dfnm]]$stations

     res[[m]][[dfnm]] <- list()
    
    # BASELINE
     
     res[[m]][[dfnm]][["baseline"]] <- estimates(ind_df, form=value~.,
                                            estimator="prequential_eval",
                                            est.pars = EST_PARS,
                                            workflow = "simple_workflow",
                                            wf.pars = WF_PARS,
                                            evaluator = "evaluate",
                                            eval.pars = EVAL_PARS,
                                            seed=SEED)

    # RANDOM RESAMPLING
     
    for(rsfun in c("under", "over")){

      for(i in 1:length(cpercs[[paste0("st", rsfun)]])){

        cperc <- cpercs[[paste0("st", rsfun)]][[i]]

        RS_PARS <- list(resample=rsfun, resample.pars = list(thr.rel=THR_REL, C.perc=cperc))
        
        parnm <- if(!is.list(cperc)) paste0(rsfun, "_cperc_", cperc) else paste0(rsfun, "_cperc_", paste0(cperc, collapse="_"))

        cat(paste("\nWith parameters", parnm))
        
         try( res[[m]][[dfnm]][[parnm]] <- estimates(ind_df, form=value~., 
                                                estimator="prequential_eval",
                                                est.pars = EST_PARS, 
                                                workflow = "simple_workflow", 
                                                wf.pars = c(WF_PARS, RS_PARS), 
                                                evaluator = "evaluate", 
                                                eval.pars = EVAL_PARS, 
                                                seed=SEED) ) 
        
      }
    }
    

    # SPATIO-TEMPORAL BIAS RESAMPLING

    # getting timestamps with right format
    if(grepl("BEIJ", dfnm)) ind_df$time <- lubridate::ymd_hms(ind_df$time)
    else ind_df$time <- as.Date(ind_df$time)
    
    for(rsfun in c("stunder", "stover")){
      
      for(i in 1:length(cpercs[[rsfun]])){
        
        cperc <- cpercs[[rsfun]][[i]]
        
        for(a in alphas){
          
          RS_PARS <- list(resample=rsfun, 
            resample.pars = list(sites_sf=stations, alpha=a, thr.rel=THR_REL, C.perc=cperc, type="add"))

          parnm <- if(!is.list(cperc)) paste0(rsfun, "_cperc_", cperc, "_alpha_", a) else paste0(rsfun, "_cperc_", paste0(cperc, collapse="_"), "_alpha_", a)
          
          cat(paste("\nWith parameters", parnm))
          
          try( res[[m]][[dfnm]][[parnm]] <- estimates(ind_df, form=value~., 
                                                 estimator="prequential_eval",
                                                 est.pars = EST_PARS, 
                                                 workflow = "simple_workflow", 
                                                 wf.pars = c(WF_PARS, RS_PARS), 
                                                 evaluator = "evaluate", 
                                                 eval.pars = EVAL_PARS, 
                                                 seed=SEED) )
          
        }
      }
    }
    
    save(res, file=paste0(RESULTS_PATH, "res_externalPrequential.Rdata"))
  }
  
}

# PROCESS AND SAVE RESULTS

library(dplyr)

sumRes <- lapply(res, function(x) lapply(x, function(y) lapply(y, function(z) z$evalRes)))

sumResTab <- bind_rows(lapply(res, function(x) 
  bind_rows(lapply(x, function(y) 
    bind_rows(lapply(y, function(z) cbind(data.frame(fold=1:nrow(z$evalRes)), 
    	as.data.frame(z$evalRes))), .id="sampling")), 
  .id="data")),
.id="model")

save(sumRes, sumResTab, file=psate0(RESULTS_PATH,"sumRes_externalPrequential.Rdata"))
