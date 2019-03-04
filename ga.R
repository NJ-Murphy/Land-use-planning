ga <- function(sparse_mat, b, constraint_ineq, obj_wts, pop_size, num_gens){
  ######################################################
  ### Code for land-use planning project
  ######################################################
  load("data/spatialDSS_inputs_fig5.RData")
  Fig5 <- as.matrix(Fig5)
  Fig5indicesfixed <- which(Fig5!=0,arr.ind = T)
  Fig5indicesother <- which(Fig5==0,arr.ind = T)
  Fig5indicesfixed <- Fig5indicesfixed[order(Fig5indicesfixed[,1]), ]
  Fig5indicesother <- Fig5indicesother[order(Fig5indicesother[,1]), ]
  
  Fig5valuesfixed <- as.vector(Fig5[Fig5indicesfixed])
  #Fig5valuesother <- as.vector(Fig5[Fig5indicesother])
  
  
  # Generate initial population of solutions (given number of items and size of popn)
  initial_pop = function(varsize,popsize){
    
    rand_vect <- function(N, M, sd = 1) {
      vec <- rnorm(N, M/N, sd)
      if (abs(sum(vec)) < 0.01) vec <- vec + 1
      vec <- round(vec / sum(vec) * M)
      deviation <- M - sum(vec)
      for (. in seq_len(abs(deviation))) {
        vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
      }
      while (any(vec < 0)) {
        negs <- vec < 0
        pos  <- vec > 0
        vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
        vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
      }
      vec
    }
    
    fixedconstraint_colinds <- numeric()
    
    for (i in 1:length(Fig5valuesfixed)){
      # the following formula gives the vector index of the grid position (index)
      fixedconstraint_colinds <- c(fixedconstraint_colinds,(Fig5indicesfixed[i,1]-1)*20*9 +  Fig5indicesfixed[i,2]*9 - (9-Fig5valuesfixed[i])) 
    }
    
    #create objective matrix
    objcoeffs <- rbind(rep(4,400),as.vector(t(Fig3a)),rep(3,400),rep(1,400),rep(5,400),rep(5,400),as.vector(t(Fig3a)),rep(7,400),as.vector(t(Fig3a)))
    Objective <- numeric()
    for (j in 1:400) {
      for (i in 1:9){
        Objective <- c(Objective,objcoeffs[i,j])
      }
    }
    #standardize each row of the matrix
    Objective <- matrix(Objective,400,9,byrow = T)
    #compute selection valie for each cell
    sigma <- t(apply(Objective,1,function(x) { (x-min(x)+0.01)/(max(x)-min(x))}))
    Objective <- numeric()
    
    pop <- matrix(0,popsize,varsize)
    pop[,fixedconstraint_colinds] <- 1
    
    for (k in 1:popsize){
      #first generate a set of possible land use types for the non-fixed cells within the min-max ranges of the no. of land types
      numvaluesforeachtype <- numeric()
      maxconstr <- numeric()
      minconstr <- numeric()
      for (i in seq(1,18,2)){
        numvaluesforeachtype <- c(numvaluesforeachtype,sample(b[i]:(b[i+1]-1),1)) #sample random number for each land type between min and max for each type
        maxconstr <- c(maxconstr,b[i+1])
        minconstr <- c(minconstr,b[i])
      }
      sizediff <- 356 - sum(numvaluesforeachtype)
      
      if (sizediff > 0){
        #add a set of possible values to add to vector numvaluesforeachtype to make sum up to 356
        numvaluesforeachtype <- numvaluesforeachtype + rand_vect(9,sizediff)
      }else if (sizediff < 0){
        #add a set of possible values to add to vector numvaluesforeachtype to make sum up to 356
        numvaluesforeachtype <- numvaluesforeachtype - rand_vect(9,-sizediff)
        
        if(any(numvaluesforeachtype<0)){ #convert ay negative values to zero
          for (i in which(numvaluesforeachtype<0)){
            sizeofnegative <- abs(numvaluesforeachtype[i])
            numvaluesforeachtype[i] <- 0
            numvaluesforeachtype[1] <- numvaluesforeachtype[1] - sizeofnegative
          }
        }
      }
      
      
      #start allocating values to each individual
      #generate matrix for fixed constraints
      fixedconstr <- rep(0,3600)
      fixedconstr[fixedconstraint_colinds] <- 1
      fixedconstr <- matrix(fixedconstr,400,9,byrow = T)
      no_landuses_fixed <- apply(fixedconstr,2,sum)
      numvaluesforeachtype <- numvaluesforeachtype + no_landuses_fixed
      
      if (sum(numvaluesforeachtype)==400){
        pre_pop <- fixedconstr
        
        #loop from 1 to 400 allocating land uses based on probabilities
        numvaluesforeachtype <- numvaluesforeachtype - no_landuses_fixed
        
        cells <- which(apply(pre_pop,1,function(x){sum(x)==0})>0) #cells to still be filled
        number_types_allocated <- rep(0,9)
        for (j in 1:356) {
          cell_allocate <- sample(cells,1)
          cells <- cells[-c(which(cells==cell_allocate))]
          
          #compute the sacling factor for each land use to encourage achievment of 
          constraintfactor <- (numvaluesforeachtype - number_types_allocated)/numvaluesforeachtype
          constraintfactor[is.nan(constraintfactor)] <- 0
          constraintfactor[is.na(constraintfactor)] <- 0
          
          value_to_allo <- sample(1:9,1,prob = sigma[cell_allocate,]*constraintfactor)
          pre_pop[cell_allocate,value_to_allo] <- 1 #allocate land use to correct column and row
          
          #print(value_to_allo)
          number_types_allocated[value_to_allo] <- number_types_allocated[value_to_allo]+1
          
        }
        if(any(apply(pre_pop,1,sum)==2)){
          toomany <- apply(pre_pop,1,sum)==2
          toolittle <- apply(pre_pop,1,sum)==0
          toomanyvec <- which(pre_pop[toomany,]>0)
          pre_pop[toomany,toomanyvec[2]] <- 0
          pre_pop[toolittle,toomanyvec[1]] <- 1
        }
        
        pop[k,] <- as.vector(t(pre_pop))
      } else {
        #generate the possible values for each of the empty cells to be allocated and then randomly scramble the values
        Fig5valuesother <- numeric()
        Fig5valuesother <- c(rep(1:9,numvaluesforeachtype))
        Fig5valuesother <- sample(Fig5valuesother,length(Fig5valuesother),replace = FALSE)
        
        othercells_colinds <- numeric()
        for (i in 1:dim(Fig5indicesother)[1]){
          othercells_colinds <- c(othercells_colinds,(Fig5indicesother[i,1]-1)*20*9 +  Fig5indicesother[i,2]*9 - (9-Fig5valuesother[i])) 
        }
        pop[k,fixedconstraint_colinds] <- 1
        pop[k,othercells_colinds] <- 1
        
      }
      
      
    }
    
    fitness_init <- evaluate_pop(pop,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds)[[1]]
    pop <- pop[-which(fitness_init<1000),] 
    
    return(list(this_pop=pop,fixedconstraint_colinds=fixedconstraint_colinds))
  }
  
  
  # Evaluate the fitness value of each member of the population 
  # Eval = 0 if no. of land types allocated constraint exceeded
  evaluate_pop = function(pop,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds){
    popsize = length(pop[,1])
    varsize = length(pop[1,])
    evals = rep(0,popsize)
    numconstraints <- dim(sparse_mat)[1]
    numconstrsatisfied <- matrix(0,popsize,1)
    min_max_constr <- matrix(0,2,9)
    min_max_constr[1,] <- b[seq(1,18,2)]
    min_max_constr[2,] <- b[seq(2,18,2)]
    
    total_obj <- apply(pop,1,function(x){sum(x*obj_wts)})
    
    for(i in 1:popsize){
      #total_constraint <- apply(sparse_mat,1,function(x){sum(pop[i,]*x)})
      
      #convert from vector to 400 by 9 matrix
      pop_temp <- matrix(pop[i,],nrow = 400,ncol = 9,byrow=TRUE)
      num_vals_per_type <- apply(pop_temp,2,sum)
      
      #check fixed constraints
      if (sum(pop[i,fixedconstraint_colinds])!=44){
        evals[i] <- 0
        #evals[i] <- total_obj[i] - 200*abs(sum(pop[i,fixedconstraint_colinds])!=44)
      } else if(sum(num_vals_per_type < min_max_constr[1,] |  num_vals_per_type > min_max_constr[2,])>0){ ## check min and max constraints
        evals[i] <- 0
        #evals[i] <- total_obj[i] - 200*abs(sum(num_vals_per_type < min_max_constr[1,] |  num_vals_per_type > min_max_constr[2,]))
      } else if (sum(pop[i,])!=400){
        evals[i] <- 0
        #evals[i] <- total_obj[i] - 200*abs(sum(pop[i,])-400)
      } else{
        evals[i] <- total_obj[i]
      }
      
    }
    return(list(evals=evals,fixedconstraint_colinds=fixedconstraint_colinds))
  }
  
  # Function for rank based selection: create 
  rank_based_selection = function(pop,fitness){
    popsize = length(pop[,1])
    varsize = length(pop[1,])
    #compute rank-based probabilities of surviving
    probs <- matrix(0,dim(pop))
    
    #get ranks of the populations fitness
    sorted_ranked_index <- order(unlist(fitness),decreasing = T, runif(length(unlist(fitness))))
    #organise population in descending order from most fit to least fit
    sorted_pop <- pop[sorted_ranked_index,]
    
    #specify size of poulation to be carried through iterations
    parent_size <- 400
    
    #linearly interpolation between 1 and xi=0.3
    probs <- approx(c(1,parent_size), y = c(1,0.5), method="linear", n=parent_size)$y
    
    #sample new population with above probabilities
    next_parents <- sorted_pop[sample(1:parent_size,size=parent_size,replace=TRUE,prob = probs),]
    
    return(next_parents)
  }
  
  # Function for crossover step (reproduction)
  reproduce = function(parents,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds){
    #children <- matrix(NA,dim(parents)[1],dim(parents)[2])
    popsize <- dim(parents)[1]
    varsize <- dim(parents)[2]
    n_rows <- 20
    n_cols <- 20
    
    #randomly generate indices of mates - mates will be successive pairs of rows
    mateind <- sample(1:popsize,popsize,replace = FALSE)
    
    # reprodution
    children <- matrix(0,dim(parents)[1],dim(parents)[2])
    for (i in seq(1,popsize,2)){
      mates <- parents[mateind[i:(i+1)],]
      
      #convert parents into 400x9 grids to find possible swaps
      mate1 <- matrix(mates[1,],nrow = 400,ncol = 9,byrow=TRUE)
      mate2 <- matrix(mates[2,],nrow = 400,ncol = 9,byrow=TRUE)
      
      mate1values <- rep(0,n_rows*n_cols)
      mate2values <- rep(0,n_rows*n_cols)
      
      #get the vlaues in each of the cells fot both parents
      mate1values <-  unlist(apply(mate1,1,function(x){as.numeric(which(x!=0))}))
      mate2values <-  unlist(apply(mate2,1,function(x){as.numeric(which(x!=0))}))
      
      land_type_combs <- combn(c(1:9),2)
      #land_type_combs <- land_type_combs[,sample(1:36,10)]
      
      for (k in 1:dim(land_type_combs)[1]){
        randswapvals <- land_type_combs[,k]
        mate2valsforindsinmate <- which((mate1values==randswapvals[1])*(mate2values==randswapvals[2])>0)
        mate1valsforindsinmate <- which((mate1values==randswapvals[2])*(mate2values==randswapvals[1])>0)
        
        #calculate hwo many of first land type in first parent and second landtype in second parent divided by 2
        numberlandtypes_for_combo_to_swap1 <- length(mate2valsforindsinmate)/2
        numberlandtypes_for_combo_to_swap2 <- length(mate1valsforindsinmate)/2
        
        #sample half of indices to be swapped
        swapinds1 <- sample(mate2valsforindsinmate,numberlandtypes_for_combo_to_swap1,replace = F)
        swapinds2 <- sample(mate1valsforindsinmate,numberlandtypes_for_combo_to_swap2,replace = F)
        
        #swap for parent 1
        mate1[swapinds1,] <- 0
        mate1[swapinds1,randswapvals[2]] <- 1
        mate1[swapinds2,] <- 0
        mate1[swapinds2,randswapvals[1]] <- 1
        
        #swap for parent 2
        mate2[swapinds1,] <- 0
        mate2[swapinds1,randswapvals[1]] <- 1
        mate2[swapinds2,] <- 0
        mate2[swapinds2,randswapvals[2]] <- 1
        
      }
      
      #convert children back to vector form
      children[i,] <- as.vector(t(mate1))
      children[i+1,] <- as.vector(t(mate2))
      
    }
    
    #childrenfit <- evaluate_pop(children,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds)
    
    ## combine parents and children and return the top 'parent size' genes from the overall population of parents and children using their fitness
    totalpop <- rbind(parents,children)
    #ranked_fitness <- as.matrix(rank(evaluate_pop(totalpop,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds)$evals,ties.method="first"))
    #sorted_ranked_index <- order(evaluate_pop(totalpop,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds)$evals,decreasing = T, runif(length(evaluate_pop(totalpop,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds)$evals)))
    
    #bestpop <- totalpop[-sorted_ranked_index[-c(1:popsize)],]
    
    return(totalpop)
  }
  
  # Function for mutation step
  mutations = function(pop,mutation_rate = 0.5){
    popsize <- dim(pop)[1]
    varsize <- dim(pop)[2]
    # MUTATE Generate the new population by mutation
    
    #which genes to mutate
    genesToMutate <- which(runif(popsize, min = 0, max = 1) < mutation_rate)
    #we then need to mutate 2 points in the gene as changing a single 1 to 0 (vice versa) will lead to violation of 1 land-use per cell
    # to do this, sample a row and cloumn index for each mutation and find another row and column in the gene with the land type you are randomly choosing to mutate the current land type to and swap these.
    pointsingenetomutaterows <- sample(1:20,length(genesToMutate),replace = T)
    pointsingenetomutatecols <- sample(1:20,length(genesToMutate),replace = T)
    
    #check if any mutations occur
    if (length(genesToMutate)!=0){
      for(i in 1:length(genesToMutate)){
        
        #convert parents into 20x20 grids of values 1 to 9 to find possible swaps
        genemutate <- matrix(pop[genesToMutate[i],],nrow = 400,ncol = 9,byrow=TRUE)
        
        genemutatevalues <- rep(0,20*20)
        for(w in 1:(20*20)){
          genemutatevalues[w] <-  as.numeric(which(genemutate[w,]!=0))
        }
        #20x20 grids of values for land-use types
        genemutategrid <- matrix(genemutatevalues,nrow = 20,ncol = 20,byrow=TRUE)
        
        #locate what value is in the row and column of the current mutation
        mutvalue <- genemutategrid[pointsingenetomutaterows[i],pointsingenetomutaterows[i]]
        #change the value to anything but its current value
        valuenew <- sample((1:9)[-mutvalue],1)
        genemutategrid[pointsingenetomutaterows[i],pointsingenetomutaterows[i]] <- valuenew
        
        #find another index in the grid that has the current land-use as this new land-use and set it to the land-use of the mutation row and column
        genemutategrid[sample(which(genemutategrid==valuenew),1)] <- mutvalue
        
        #convert grids back into vectors of 3600 but first set the current gene thats being mutated to zeros
        pop[i,] <- 0
        for (m in 1:20){
          for (n in 1:20){
            pop[i,(m-1)*20*9 + (n-1)*9 + genemutategrid[m,n]] <- 1
          }
        }
      }
    }
    return(pop)
  }
  
  ## putting it all together
  varsize <- length(obj_wts)
  this_pop <- initial_pop(varsize,pop_size)
  fixedconstraint_colinds <- this_pop$fixedconstraint_colinds
  this_pop <- this_pop$this_pop
  
  max_evals = c()
  mean_evals = c()
  nconstraintssat = c()
  
  for(gen in 1:num_gens){
    evals = evaluate_pop(this_pop,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds)
    next_parents = rank_based_selection(pop=this_pop,fitness=evals$evals)
    next_children = reproduce(parents=next_parents,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds)
    next_children = rank_based_selection(pop=next_children,fitness=evaluate_pop(next_children,obj_wts,sparse_mat,b,constraint_ineq,fixedconstraint_colinds)$evals)
    this_pop = mutations(pop=next_children,mutation_rate=0.05) 
    #this_pop = next_children
    max_evals = c(max_evals,max(evals$evals))
    mean_evals = c(mean_evals,mean(evals$evals))
    #nconstraintssat <- c(nconstraintssat,max(evals[[2]]))
  }
  
  
  return(list(evals=evals$evals,pop=this_pop,maxeval=max_evals,meaneval=mean_evals,ngenerations = num_gens))
}
## EOF - GA