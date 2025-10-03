################################################################################
## SOURCE: https://github.com/gtonkinhill/panstripe/blob/06bc5e9cca7e9f9eb74c89a8cb5247d6983246f2/R/asr_parsimony.R

# This code is a modified version of a function from the castor package originally written by Stilianos Louca (https://cran.r-project.org/web/packages/castor/)
# Maximum parsimony ancestral state reconstruction for discrete traits.
# Modification of Sankoff algorithm for reconstructing discrete ancestral states (Weighted Small Parsimony Problem)
# Sankoff's algorithm allows the inclusion of a cost matrix:
#  	transition_costs[i,j] is the cost of transitioning i-->j (ignoring edge length)
# 	If transition_costs is "all_equal", then all transitions are penalized equally (same as if transition_costs[i,j] = 1-delta_{ij})
# 	If transition_costs is "sequential", then only single-step transitions (i-->i+1) are allowed, and all are penalized equally
# 	If transition_costs is "proportional", then all transition are allowed but they are penalized proportionally to the number of steps.
#  The modification of this function is that optionally, edge lengths can be used to weight the transition costs:
#  	Longer edges imply smaller transition costs between states
#  	Specifically, the cost of transitioning is transition_cost[i,j]/(edge_length^edge_exponent)
# 	where edge_exponent can be e.g. 0 (Sankoff's original algorithm), 1 (linear weighting) or 0.5 (square-root weighting, corresponding to a Brownian motion)
#  Requirements:
# 	Tree can be multifurcating, and can also include nodes with a single child
#	Tree must be rooted.
# 	If (edge_exponent>0) then: All edges must have length > 0
#  For a description of the original Sankoff algorithm, see: 
# 	http://telliott99.blogspot.ca/2010/03/fitch-and-sankoff-algorithms-for.html
# 	(page 11) https://cs.brown.edu/courses/csci1950-z/slides/CSCI1950ZFall09_Lecture2.pdf
#  The function returns ancestral state probabilities as a (non-flattened) NumericMatrix of size Nnodes x Nstates.
asr_max_parsimony = function(	tree, 
                              tip_states, 			# integer vector of size Ntips
                              Nstates				= NULL, 
                              transition_costs	= "all_equal", 
                              edge_exponent		= 0,
                              weight_by_scenarios	= TRUE,
                              check_input			= TRUE){
  Ntips  = length(tree$tip.label)
  Nedges = nrow(tree$edge)
  
  # basic error checking
  if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
  if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be integers"))	
  if(is.null(Nstates)) Nstates = max(tip_states);
  if(check_input){
    min_tip_state = min(tip_states)
    max_tip_state = max(tip_states)
    if((min_tip_state<1) || (max_tip_state>Nstates)) stop(sprintf("ERROR: tip_states must be integers between 1 and %d, but found values between %d and %d",Nstates,min_tip_state,max_tip_state))
    if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
  }
  
  # construct transition_costs matrix if needed
  if(is.character(transition_costs)){
    if(transition_costs=="all_equal"){
      # all transitions penalized equally
      transition_costs = matrix(1, nrow=Nstates, ncol=Nstates)
    }else if(transition_costs=="sequential"){
      # only single-step transitions are allowed, and all are penalized equally
      transition_costs = matrix(Inf, nrow=Nstates, ncol=Nstates)
      for(i in 1:Nstates){
        if(i<Nstates) transition_costs[i,i+1] = 1
        if(i>1) transition_costs[i,i-1] 	  = 1
      }
    }else if(transition_costs=="proportional"){
      # all transitions are allowed, but penalized proportional to the number of steps
      transition_costs = matrix(0, nrow=Nstates, ncol=Nstates)
      for(i in 1:Nstates){
        transition_costs[i,] = sapply(1:Nstates, function(j) abs(j-i))
      }
    }else if(transition_costs=="exponential"){
      # all transitions are allowed, but penalized exponentially to the number of steps
      transition_costs = matrix(0, nrow=Nstates, ncol=Nstates)
      for(i in 1:Nstates){
        transition_costs[i,] = sapply(1:Nstates, function(j) exp(abs(j-i)))
      }
    }else{
      stop(sprintf("ERROR: Uknown transition_costs '%s'",transition_costs));
    }
    diag(transition_costs) = 0; # no cost for staying in the same state
  }else{
    if(nrow(transition_costs)!=Nstates || ncol(transition_costs)!=Nstates) stop(sprintf("ERROR: Transition costs has wrong size (%d x %d), expected %d x %d",nrow(transition_costs),ncol(transition_costs),Nstates,Nstates));
    if(check_input){
      if(any(transition_costs<0)) stop(sprintf("ERROR: Some transition costs are negative (found value %g)",min(transition_costs)))
      if(any(diag(transition_costs)!=0)) stop(sprintf("ERROR: The diagonal of transition_costs includes non-zero values, which makes no sense in a maximum-parsimony model"))
      if((edge_exponent!=0) && (!is.null(tree$edge.length)) && (any(tree$edge.length==0))) stop(sprintf("ERROR: edge_exponent is non-zero, but some edges in the tree have zero length"))
    }
  }
  
  ##############################################################################
  ## Start edit
  
  ## add panstripe:::
  results = panstripe:::WMPR_ASR_CPP(	Ntips			 						= Ntips,
                                      Nnodes			 						= tree$Nnode,
                                      Nedges			 						= Nedges,
                                      Nstates			 						= Nstates,
                                      tree_edge 		 						= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
                                      edge_length		 						= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
                                      tip_states		 						= tip_states-1,
                                      transition_costs 						= as.vector(t(transition_costs)),
                                      branch_length_exponent 					= edge_exponent,
                                      weight_posteriors_by_scenario_counts	= weight_by_scenarios,	# (INPUT) if true, then the posterior_probability of a state (in a specific node) is proportional to the number of scenarios in which the trait is at that state
                                      verbose									= FALSE,
                                      verbose_prefix							= "");
  
  # 重构 C++ 输出
  scenario_counts   <- matrix(results$scenario_counts,   ncol=Nstates, byrow=TRUE)  # Nnodes × Nstates
  transition_counts <- matrix(results$transition_counts, byrow=TRUE, ncol=2)        # Nnodes × 2，col1=0→1，col2=1→0
  
  # 找到每个 tip 的父节点在 internal-node 序号中的行号
  is_above_tip <- tree$edge[,1][ match(1:Ntips, tree$edge[,2]) ]
  ia <- is_above_tip - Ntips   # 转为 internal-node 索引（1…Nnodes）
  
  # 总共的“边”数：tips + internals
  N_total <- Ntips + tree$Nnode
  
  # 初始化 gain/loss 向量
  gain <- numeric(N_total)
  loss <- numeric(N_total)
  
  # （1）内部节点：期望事件数 = 各类型事件计数 / 场景总数
  #    transition_counts[,1] = 0→1 次数, transition_counts[,2] = 1→0 次数
  sc_sum <- rowSums(scenario_counts)   # 每个 internal node 的场景数总和
  internals <- (Ntips+1):N_total
  gain[internals] <- transition_counts[,1] / sc_sum
  loss[internals] <- transition_counts[,2] / sc_sum
  
  # （2）叶节点：把原来 changes 的逻辑拆成 gain vs loss
  #    posterior_probabilities[row,1] = P(parent=0)
  #    posterior_probabilities[row,2] = P(parent=1)
  #  tip_states 原始取值为 1 或 2，对应 child state 0/1
  gain[1:Ntips] <- results$posterior_probabilities[ia, 1] * (tip_states == 2)
  loss[1:Ntips] <- results$posterior_probabilities[ia, 2] * (tip_states == 1)
  
  # 最后 return 时分别输出
  return(list(
    success               = TRUE,
    ancestral_likelihoods = results$posterior_probabilities,
    scenario_counts       = scenario_counts,
    total_cost            = results$best_root_cost,
    changes_gain          = gain,   # 每条边上 0→1（基因获得）的期望次数/概率
    changes_loss          = loss    # 每条边上 1→0（基因丢失）的期望次数/概率
  ));
}


################################################################################
## https://github.com/gtonkinhill/panstripe/blob/06bc5e9cca7e9f9eb74c89a8cb5247d6983246f2/R/panstripe.R

twd_llk <- function(p, model, data) {
  # suppressWarnings({
  tm <- stats::glm(model, data, family = statmod::tweedie(var.power = p, link.power = 0))
  # })
  tsm <- summary(tm)
  -sum(log(tweedie::dtweedie(y = tm$y, mu = tm$fitted.values, phi = tsm$dispersion, power = p)))
}

fit_tweedie <- function(model, data, method='glm'){
  
  if (method=='glm'){
    fm <- tryCatch(
      {
        op <- stats::optimise(twd_llk, lower = 1, upper = 2, model=model, data=data)
        tm <- stats::glm(model, data, family = statmod::tweedie(var.power = op$minimum, link.power = 0))
        stm <- summary(tm)
        tm$p <- op$minimum
        tm$phi <- stm$dispersion
        tm
      },
      error=function(cond) {
        stop(
          "Panstripe model fit failed! This can sometime be caused by unusual branch lengths.
Setting fit_method='glmmTMB' or family='quasipoisson' or 'gaussian' often provides a more stable fit to difficult datasets"
        )
      }
    )
  } else {
    fm <- glmmTMB::glmmTMB(model, data = data, family = glmmTMB::tweedie)
    fm$coefficients <- fm$fit$par[1:(length(fm$fit$par)-2)]
    fm$p <- min(1.99, max(1.01, exp(1+fm$fit$par[[length(fm$fit$par)]])))
    fm$phi <- exp(fm$fit$par[[length(fm$fit$par)-1]])
    
    if (fm$fit$convergence!=0){
      warning(
        "Panstripe model fit failed to converge!
Setting family='quasipoisson' or 'gaussian' often provides a more stable fit to difficult datasets"
      )
    }
  }
  
  return(fm)
}

fit_model <- function(d, indices=1:nrow(d), tree=NULL, min_depth=NULL, model, family, fit_method, boot_type){
  stopifnot(length(indices)==nrow(d))
  stopifnot(boot_type %in% c('branch', 'gene'))
  
  # sometimes it is hard to get the  model to converge for a particular sample. We attempt a few which hopefully does not bias things too much.
  coef <- NULL
  attempt <- 0
  max_attempt <- 5
  
  while(is.null(coef) && (attempt<=max_attempt)){
    attempt <- attempt + 1
    try({
      if (boot_type == 'branch') {
        tdat <- d[indices,]
      } else {
        tdat <- tibble::tibble(
          acc=colSums(d[indices,]),
          core=tree$edge.length,
          istip=tree$edge[,2]<=length(tree$tip.label) 
        )
        # add depth
        tdat$depth <- ape::node.depth.edgelength(tree)[tree$edge[,2]]
        if (!is.null(min_depth)){
          tdat <- tdat[tdat$min_depth>=min_depth,] 
        }
      }
      
      if (is.character(family) && (family=="Tweedie")){
        tm <- fit_tweedie(model, data = tdat, method=fit_method)
        coef <- c(tm$coefficients, tm$p, tm$phi)
      } else {
        if (fit_method=='glm'){
          tm <- stats::glm(model, data = tdat, family=family)
        } else {
          tm <- glmmTMB::glmmTMB(model, data = tdat, family=family)
          tm$coefficients <- tm$fit$par[names(tm$fit$par)=='beta']
        }
        coef <- c(tm$coefficients, NA, NA)
      }
    })
  }
  
  if (is.null(coef) && (attempt==max_attempt+1)){
    stop('Model fitting failed to converge in bootstrap replicate!')
  }
  return(coef)
}

boot_ci_pval <- function(boot_res, index, type, 
                         theta_null=0, 
                         precision=NULL, 
                         ci_conf=0.95,
                         transformation='identity',
                         calc_pval=FALSE) {
  
  if (!transformation %in% c('identity','logit','log')) stop('Invalid transformation!')
  
  if (is.null(precision)){
    precision = 1/boot_res$R  
  }
  
  if (transformation=='logit'){
    ll <- stats::make.link('logit')
    h <- function(x) ll$linkfun(x-1)
    hinv <- function(x) ll$linkinv(x)+1
    hdot <- function(x) ll$mu.eta(x)
  } else if (transformation=='log'){
    ll <- stats::make.link('log')
    h <- function(x) ll$linkfun(x)
    hinv <- function(x) ll$linkinv(x)
    hdot <- function(x) ll$mu.eta(x)
  } else {
    h <- function(x) x
    hinv <- function(x) x
    hdot <- function(x) 1
  }
  
  if (calc_pval){
    #calculate p-value by inverting the corresponding confidence intervals, as described in Section 3.12 in Hall (1992)
    alpha_seq <- seq(precision, 1 - 1e-16, precision)
    ci <- suppressWarnings({boot::boot.ci(boot_res, conf = 1 - alpha_seq, 
                                          type = type, index=index,
                                          h = h, hinv = hinv, hdot = hdot)})
    bounds <- switch(type, 
                     norm = ci$normal[, 2:3],
                     basic = ci$basic[,4:5],
                     stud = ci$student[, 4:5],
                     perc = ci$percent[, 4:5], 
                     bca = ci$bca[, 4:5])
    alpha <- alpha_seq[which.min(theta_null >= bounds[, 1] & 
                                   theta_null <= bounds[, 2])]
  }
  
  #calculate CI
  ci <- suppressWarnings({boot::boot.ci(boot_res, conf = ci_conf, 
                                        type = type, index=index,
                                        h = h, hinv = hinv, hdot = hdot)})
  bounds <- switch(type, 
                   norm = ci$normal[, 2:3],
                   basic = ci$basic[,4:5],
                   stud = ci$student[, 4:5],
                   perc = ci$percent[, 4:5], 
                   bca = ci$bca[, 4:5])
  
  if (calc_pval){
    return(c(sort(bounds), alpha))
  } else {
    return(sort(bounds))
  }
}
