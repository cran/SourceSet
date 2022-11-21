#' Simulated dataset
#'
#' This data contains the parameters used in the study of the finite case behavior of source set algorithm, as described in of Salviato et al. (2019).
#' @docType  data
#' @usage data(simulation)
#' @format A list class that contains the true parameters (\code{mu}, vector of means and \code{S}, covariances matrix) of two multivariate normal distributions in two different experimental conditions (\code{condition1}, reference condition and \code{condition2}, perturbed condition) and the underlying graphical structure G (\code{graph}. Six different perturbations are considered, see below. ).
#'
#' The differences between the two conditions are driven by:
#' \itemize{
#'   \item{ a node that is a separator within the graph (\code{simulation$condition2$`5`})}
#'   \item{ a node that is contained in only one clique of the graph (\code{simulation$condition2$`10`})}
#' }
#'
#' The intensity of the artificial perturbation is:
#' \itemize{
#'   \item{ mild (\code{simulation$condition2$`10`$`1.2`})}
#'   \item{ moderate (\code{simulation$condition2$`10`$`1.6`})}
#'   \item{ strong (\code{simulation$condition2$`10`$`2`})}
#' }
#' @details The starting parameters of the reference condition are obtained by randomly selecting a gene set of the same cardinality as the order of the graph G, from the Acute Lymphocytic Leukemia (\code{\link[ALL]{ALL}}) dataset. These are then modified to represent the parameters of the perturbed condition.
#' Formally, starting from the parameters related to the reference group, the procedure act on means and variances so that the conditional distribution of the variables on which it does not directly intervene remains unchanged under the two conditions.
#' However, this action affects the entire global joint distribution, thus creating the propagation effect. See Salviato et al. (2016) for more details.
#' @seealso \code{simPATHy}, \code{\link[ALL]{ALL}}
#' @references
#' Chiaretti, S. et al. (2005). Gene expression profiles of b-lineage adult acute lymphocytic leukemia reveal genetic patterns that identify lineage derivation and distinct mechanisms of transformation. Clinical Cancer Research, 11(20), 7209–7219.
#'
#' Salviato, E. et al. (2016). \code{simPATHy}: a new method for simulating data from perturbed biological pathways. Bioinformatics, 33(3), 456–457.
#'
#' Salviato, E. et al. (2019). \code{SourceSet}: a graphical model approach to identify primary genes in perturbed biological pathways. PLoS computational biology 15 (10), e1007357.
"simulation"


#' Default shrinkage estimation of covariance matrices
#'
#' The function adds a small quantity to the diagonals of covariance matrix estimates to regularize them.
#' @param s covariance matrix estimte in the pooled sample
#' @param s1 sample covariance matrix estimate in class 1
#' @param s2 sample covariance matrix estimate in class 2
#' @param param list of parameters: \itemize{\item{\code{lamda}: a vector of lambdas (to be supplied only if some custom lambdas are to be used)} \item{\code{type}: minimum (\code{min}), maximum (\code{max}) or optimal (\code{opt})} \item {\code{probs}: the numeric value of probability with value in [0,1]}}
#' @details To determine the quantity to add to the diagonals of covariance matrices, the function:
#' \itemize{
#'  \item{ finds the distributions of the sample variances of the p variables in the two classes and in the pooled sample}
#'  \item{ computes the \code{probs} percentile of each of these distributions }
#'  \item{ use the minimum, maximum or optimal (one for each matrix) (\code{type})}
#' }
#' @note It should be stressed that the default parameters for TEGS shrink estimator allow to compare the log likelihood criterion among distributions if the \code{\link[SourceSet]{testMeanVariance}} is performed.
#' @references Huang, Y.-T. and Lin, X. (2013). Gene set analysis using variance component tests. BMC Bioinformatics, 14(1), 210.
#' @seealso \code{\link[SourceSet]{testMeanVariance}}, \code{\link[SourceSet]{parameters}}
#' @return The function returns a list of shrink covariance matrices.
#' @examples
#'
#'
#' if(require(mvtnorm)){
#'
#'   ## Generate two random samples of size 50 from two multivariate normal distributions
#'   # sample size
#'   n<-50
#'   # true parameters of class 1 and class 2
#'   param.class1<-simulation$condition1
#'   param.class2<-simulation$condition2$`5`$`2`
#'   # simulated dataset
#'   data.class1<-rmvnorm(n = n,mean =param.class1$mu ,sigma =param.class1$S)
#'   data.class2<-rmvnorm(n = n,mean =param.class2$mu ,sigma=param.class2$S)
#'   data<-rbind(data.class1,data.class2)
#'   classes<-c(rep(1,nrow(data.class1)),rep(2,nrow(data.class2)))
#'
#'   ## estimated parameters: maximum likelihood estimate
#'   s<-cov(data)
#'   s1<-cov(data.class1)
#'   s2<-cov(data.class2)
#'
#'   ## default parameters:
#'   # use the minimum of median variances distributions of the three supplied covariance matrices
#'   def.shrink<-shrinkTEGS(s,s1,s2)
#'   def.shrink$lambda
#'
#'   ## use customize lamdas
#'   def.shrink<-shrinkTEGS(s,s1,s2,param = list(lambda=c(0.1,0.2,0.3)))
#'   def.shrink$lambda
#'
#'   # use for each covariance matrix the 0.4 percentile of its variances distributions
#'   def.shrink<-shrinkTEGS(s,s1,s2,param = list(type="opt",probs=0.4))
#'   def.shrink$lambda
#' }
#'
#' @export
shrinkTEGS<-function(s,s1,s2,param=list(probs=0.05,type="min")){

  ## check lambda
  if(!is.null(param$lambda) & length(param$lambda)==3){
    lambda<-param$lambda
    param$type<-"user"

  } else {
    fun.type<-switch(param$type,
                     min=function(x){ return(rep(min(x),3)) },
                     max=function(x){ return(rep(max(x),3)) },
                     opt=function(x){ return(x) }
    )
    lambda<-fun.type(sapply(list(diag(s),diag(s1),diag(s2)),stats::quantile,probs=param$probs))
  }

  p<-ncol(s)

  # add to the diagonal the variance quantile
  res<-list(S=s+diag(p)*lambda[1],
            S1=s1+diag(p)*lambda[2],
            S2=s2+diag(p)*lambda[3],
            lambda=list(lambda=lambda,type=param$type,probs=param$probs)
            )

  return(res)
}

#' Estimation of parameters for test equality of two normal distributions
#'
#' The function estimates the parameters of two normal distributions. Both maximum likelihood estimates and shrinkage estimate of covariance matrices are supplied.
#' @param data an expression matrix with colnames for variables and row names for samples
#' @param classes a vector of length equal to the number of rows of data. It indicates the class (condition) of each statistical unit. Only two classes, labeled as 1 and 2, are allowed
#' @param shrink boolean. if \code{FALSE} the maximum likelihood estimates are returned; if \code{TRUE} the shrinkage estimates are returned instead
#' @param shrink.function function that implements the shrinkage method. It must return a list object with all the elements required as input arguments in \code{\link[SourceSet]{testMeanVariance}}. Default is \code{\link[SourceSet]{shrinkTEGS}} function.
#' @param shrink.param additional parameters to pass as input arguments of the shrink function specified in the \code{shrink.function}.
#' @seealso \code{\link[SourceSet]{shrinkTEGS}}, \code{\link[SourceSet]{testMeanVariance}}
#' @return The function returns a list containing: three matrices with maximum likelihood or the shrinkage estimates (pooled, condition1 and condition2), and a list with the used parameters.
#' @examples
#' if(require(mvtnorm)){
#'   ## Generate two random samples of size 50 from two multivariate normal distributions
#'   # sample size
#'   n<-50
#'   # true parameters of class 1 and class 2
#'   param.class1<-simulation$condition1
#'   param.class2<-simulation$condition2$`5`$`2`
#'   # simulated dataset
#'   data.class1<-rmvnorm(n = n,mean =param.class1$mu ,sigma =param.class1$S)
#'   data.class2<-rmvnorm(n = n,mean =param.class2$mu ,sigma=param.class2$S)
#'   data<-rbind(data.class1,data.class2)
#'   classes<-c(rep(1,nrow(data.class1)),rep(2,nrow(data.class2)))
#'
#'   ## estimated parameters: maximum likelihood estimate
#'   est.param<-parameters(data = data,classes =classes ,shrink = FALSE)
#'
#'   ## estimated parameters: regularized estimate
#'   est.param.shrink<-parameters(data = data,classes =classes ,shrink = TRUE)
#'   # tuning values and other info on shrinkage estimate
#'   str(est.param.shrink$shrink.info)
#' }
#' @export
parameters<-function(data, classes, shrink =TRUE,shrink.function=shrinkTEGS,shrink.param=list(probs=0.05,type="min")){

  n1<-sum(classes==1)
  n2<-sum(classes==2)

  d1<-data[classes==1,,drop=FALSE]
  d2<-data[classes==2,,drop=FALSE]

  p<-ncol(data)

  s<-stats::cov(data)
  s1<-stats::cov(d1)
  s2<-stats::cov(d2)


  if(shrink) estim<-shrink.function(s = s,s1 = s1,s2 = s2,param = shrink.param)
  else estim<-list(S=s,S1=s1,S2=s2,lambda=NULL)


  res<-list(
    S=estim$S,S1=estim$S1,S2=estim$S2,
    n1=n1,n2=n2,
    shrink.info=list(shrink=shrink,lambda=estim$lambda)
  )

  return(res)

}

#' Test the equality of two normal distributions
#'
#' The function performs the test of equality of two multivariate normal distrbutions (class1 and class2).
#' @param S estimated covariance matrix for pooled sample
#' @param S1 estimated covariance matrix in class 1
#' @param S2 estimated covariance matrix in class 2
#' @param n1 number of samples in class 1
#' @param n2 number of samples in class 2
#' @details The criterion for testing the equality of two normal distributions is the following:
#' \deqn{\Lambda_c =  n_1 * log( |S| / |S^1| ) +  n_2 * log( |S| / |S^2| )  }
#' The asymptotic null distribution of the criterion, when the maximum likelihood estimates of the covariance matrices are used, is Chi square with \eqn{ |\Gamma| * (|\Gamma|+3) / 2} degrees of freedom, where G is the dimension of the underlying distributions.
#' @note The asymptotic null distributions holds only when the maximum likelihood estimates of the covariance matrices are supplied.
#' @seealso \code{\link[SourceSet]{parameters}}
#' @examples
#' if(require(mvtnorm)){
#'
#'   ## Generate two random samples of size 50 from two multivariate normal distributions
#'   # sample size
#'   n<-50
#'   # true parameters of class 1 and class 2
#'   param.class1<-simulation$condition1
#'   param.class2<-simulation$condition2$`5`$`2`
#'   # simulated dataset
#'   data.class1<-rmvnorm(n = n,mean =param.class1$mu ,sigma =param.class1$S)
#'   data.class2<-rmvnorm(n = n,mean =param.class2$mu ,sigma=param.class2$S)
#'   data<-rbind(data.class1,data.class2)
#'   classes<-c(rep(1,nrow(data.class1)),rep(2,nrow(data.class2)))
#'
#'   s<-cov(data)
#'   s1<-cov(data.class1)
#'   s2<-cov(data.class2)
#'   testMeanVariance(S = s,S1 =s1, S2 = s2, n1 = n, n2 = n)
#'
#'   ## equivalently...
#'   # estimated parameters: maximum likelihood estimate
#'   est.param<-parameters(data = data,classes =classes ,shrink = FALSE)
#'   testMeanVariance(est.param$S,est.param$S1,est.param$S2,est.param$n1,est.param$n2)
#' }
#' @return The function returns a list that contain the test statistic (\code{stat}) and the p-value test obtained of equality, using the asymptotic distribution (\code{alpha}).
#' @export
testMeanVariance<-function(S,S1,S2,n1,n2){

  p<-dim(S)[1]
  gdl<- p+(p*(p+1))/2


  s1D <- det( S1* ((n1-1)/n1) )
  s2D <- det( S2* ((n2-1)/n2 ) )
  sD<- det( S* ((n1+n2-1)/(n1+n2)) )

  lc <- n1 * log(sD/s1D) + n2 * log(sD/s2D)
  ac<-  1-stats::pchisq(lc,gdl )

  return(list(stat=lc,alpha=ac))
}



#' All possible RIP orderings
#'
#' The function identifies all possible clique orderings leading to distinct factorizations of the associated joint distribution.
#' @param graph a graph represented as a \code{graphNEL} object. If the input graph is not decomposable, the function will internally moralize and triangulate it.
#' @details For each root clique, the function uses \code{\link[gRbase]{rip}} function to identify a sequence of the set of cliques that satisfies the running intersection property by first ordering variables by the maximum cardinality search algortithm. The \code{root} argument is used to check which clique will be the first to enter in the rip ordering.
#' @seealso \code{\link[gRbase]{rip}}
#' @return Given a graph, the function returns:
#' \itemize{
#'  \item{\code{elements}:  a list composed of four other lists:
#'   \itemize{
#'     \item{\code{cliques}: (a list of character vectors) variables contained in each maximal clique of the moralized and triangulated input graph}
#'     \item{\code{separators}: (a list of character vectors) unique separators, i.e., common variables among cliques}
#'     \item{\code{components}: (a dataframe) unique “clique | separator” elements defined on the basis of all rip orderings. Each element represents a conditional distribution (see Djordjilovic and Chiogna)}
#'     \item{\code{variables}: (a character vector) nodes of the graph }
#'   }
#'  }
#'  \item{\code{indices}: a list composed of two other lists:
#'    \itemize{
#'     \item{ \code{all}: (a list of character vectors) \code{cliques} and \code{separators} lists of variables}
#'     \item{ \code{ordering}: one dataframe for each identified ordering. Each data frame is a subset of size k (i.e., number of maximal cliques), of the \code{components} elements. The name of each list corresponds to the used root clique.}
#'    }
#'  }
#'  \item{\code{graph}: decomposable graph used in the identification of rip orderings. It may differ from the input graph. In fact, if the input graph is not decomposable, the function will internally moralize and triangulate it.}
#' }
#' @examples
#' if(require(gRbase)){
#'
#'   ## decomposable graph
#'   ug.graph<-ug(~1:2:3+3:4+4:5:6:7)
#'   ug.rip.all<-ripAllRootsClique(ug.graph)
#'   # 7 variables
#'   length(ug.rip.all$elements$variables)
#'   # 3 max.cliques
#'   length(ug.rip.all$elements$cliques)
#'   # 7 unique components
#'   nrow(ug.rip.all$elements$components)
#'   # all rip orderings:
#'   ug.rip.all$indices$ordering
#'
#'   ## directed graph
#'   dag.graph<-dag(~3:1+3:2+4:3)
#'   dag.rip.all<-ripAllRootsClique(dag.graph)
#'   # triangulated and morliazed graph
#'   dag.rip.all$graph
#'   # all rip orderings
#'   dag.rip.all$indices$ordering
#' }
#' @export
ripAllRootsClique<-function(graph){


  ## check graph
  if( !inherits(graph,"graphNEL") ){ #!(class(graph)=="graphNEL")
     stop("graph argument is not a graphNEL object")
   } else {

     ## Internal function
     graph<-moralize(graph)
     #if(gRbase::is.DG.graphNEL(graph) ) {
     #  #graph<- gRbase::moralize(graph)
     #   graph<- clipper:::mmmoralize(graph)
     #  message(".. moralize graph\n")
     #} else { graph<- clipper:::mmmoralize(graph) }

     if( !gRbase::is.TUG(graph) ) {
       graph<- gRbase::triangulate(graph)
       #message(".. triangulate graph\n")
     }
   }

  variables<-graph::nodes(graph)

  # internal function to convert names
  convert.name<-function(name,lconv){
    ind<-match(lconv,name)
    names(name)[ind]
  }

  ## get cliques (must be unique)
  cliques<-gRbase::getCliques(graph)
  cliques<-lapply(cliques,sort)
  names(cliques)<-paste0("C",1:length(cliques))

  ## set each clique as root
  ripAll<-lapply(cliques,function(cl,gr) gRbase::rip(object = gr,root =cl )  ,gr=graph)

  ## get separators (could be not unique)
  # the empty separator is always S0
  separators<-unique(ripAll[[1]]$separators)
  separators<-lapply(separators,sort)
  names(separators)<-paste0("S",0:c(length(separators)-1))

  ## convert ripAll names: they must match with orginal cliques
  # match works if the elements are identical (sort all the names)
  for(i in 1:length(ripAll)){
    # clique name conversion
    ripAll[[i]]$cliques<-lapply(ripAll[[i]]$cliques,sort)
    names(ripAll[[i]]$cliques)<-convert.name(cliques,ripAll[[i]]$cliques)

    # separator name conversion
    ripAll[[i]]$separators<-lapply(ripAll[[i]]$separators,sort)
    names(ripAll[[i]]$separators)<-convert.name(separators,ripAll[[i]]$separators)
  }

  ## find components
  # the last columns is S0
  ind<-which(names(separators)=="S0")
  all<-append(cliques,separators[-ind])
  all<-append(all,separators[ind])


  ordering<-lapply(ripAll,function(x,n){
    df<-data.frame(clique=names(x$cliques),
                   separator=names(x$separators),
                   stringsAsFactors = FALSE)

    df<-data.frame(df,ind_clique=match(df$clique,n),
                   ind_separator=match(df$separator,n))

  },n=names(all))

  components<-unique(Reduce(rbind,ordering))
  components<-data.frame(components,component=paste0("comp",1:nrow(components)),stringsAsFactors = FALSE)

  ordering<-lapply(ordering,function(x,y) merge(x,y),y=components)

  elements<-list(cliques=cliques,separators=separators,components=components,variables=variables)
  indices<-list(all=all, ordering=ordering)

  ordering<-list(elements=elements,indices=indices,graph=graph)
  class(ordering)<-"allOrd"

  return( ordering )
}

#' Get random permutations of a set of elements
#'
#' The function arranges, in an optimized way, all the elements of a set into a selected number of different sequences (i.e., permutations). If the number of possible orderings is less than the required number, the function returns the collection of all possible permutations.
#' @param n number of elements
#' @param nperms number of required permutations
#' @return
#' The function returns:
#' \itemize{
#'    \item{ \code{perms}: a matrix with \code{nperms} rows and \code{n} columns, containing the sequence of the ordered elements}
#'    \item{ \code{all.perms.flag}: \code{1} if the perms array contains the entire collection of permutations, \code{0} otherwise. In the first case, the number of rows of \code{perms} matrix may be less than the number of requested permutations}
#'    \item{ \code{nperms.act}: the number of permutations returned}
#' }
#' @examples
#' sub.perm<-getPermutations(10,100)
#' all.perm<-getPermutations(3,100)
#' @export
getPermutations<-function (n, nperms) {
  if (n > 50)
    total.perms <- Inf
  else total.perms = factorial(n)
  if (total.perms <= nperms) {
    perms = gtools::permute(1:n)
    all.perms.flag = 1
    nperms.act = total.perms
  }
  if (total.perms > nperms) {
    perms = matrix(NA, nrow = nperms, ncol = n)
    for (i in 1:nperms) {
      perms[i, ] = sample(1:n, size = n)
    }
    all.perms.flag = 0
    nperms.act = nperms
  }
  return(list(perms = perms, all.perms.flag = all.perms.flag,
              nperms.act = nperms.act))
}

# Get the Mean Variance test distribution for each clique
#
# Calculate the asympotic or the permutational distribution of the MeanVariance test
# @param data ...
# @param classes ...
# @param permute ...
# @param nperms ...
# @param indcs ...
# @param shrink ...
# @param seed ...
# @param shrink.function ...
# @param shrink.param ...
# @param print.label ...
testMeanVarianceDistribution<-function(data,classes,
                                       permute=TRUE,nperms=1000,
                                       indcs,shrink=FALSE,seed=NULL,
                                       print.label=NULL,shrink.function=shrinkTEGS, shrink.param=list(probs=0.05,type="min")){

  n<-nrow(data)
  p<-max(sapply(indcs,length))
  nc<- min(table(classes))

  ## check number of classes
  if( length(table(classes))!=2 ){
    stop("the number of distinct classes must to be two")
  }

  ## check permute and shrink
  if( shrink & !permute){
    warning("permute argument is not valid. The asymptotic distribution can't be used for the shrinkage estimate of the covariance matrices. permute=TRUE option is used.")
    permute<-TRUE
  }

  ## Check if it is possible to use the asympotic test
  # Case: shrink=FALSE and nc<=p
  if( shrink == FALSE  & nc<= p ){
    msg<-paste0("The asymptotic mean-variance test can't be used. The number of samples in the smallest class (i.e. ",nc,") must to be greater then number of variables in the maximal clique (i.e. ",p,"): shrink=TRUE option is used.")
    warning(msg)
    shrink<-TRUE
  }

  if(!is.null(seed)) set.seed(seed)
  #index<-lapply(elem,function(x,n){ match(x,n) },n=colnames(data)  )
  index<-indcs

  if(permute) mat<-rbind(1:nrow(data),getPermutations(n,nperms)$perms)
  else mat<-matrix(1:nrow(data),nrow = 1)

  matDist<-matrix(NA,ncol=length(index),nrow=nrow(mat))
  colnames(matDist)<-names(index)
  rownames(matDist)<-0:c(nrow(mat)-1)

  N<-nrow(mat)

  ## compute unique lambda
  ObsParam<-parameters(data = data,classes = classes,shrink = shrink, shrink.function = shrink.function, shrink.param = shrink.param )
  shrink.param_i<-ObsParam$shrink.info$lambda

  if(permute) {
    msg<-paste0("( ",print.label," )"," permutations [:bar] :percent in :elapsed")
    pb<- progress::progress_bar$new(
      format=msg,
      total=N,clear=FALSE,width=80
    )
  }

  for(i in 1:N){

    data_i<-data[mat[i,],,drop=FALSE]
    rownames(data_i)<-cl_i<-paste(classes)

    param_i<-parameters(data = data_i,classes = cl_i,shrink = shrink, shrink.function = shrink.function, shrink.param = shrink.param_i )

    dist_i<-sapply(index,function(x,p){
      testMeanVariance(S = p$S[x,x,drop=FALSE],S1 =p$S1[x,x,drop=FALSE] ,S2 = p$S2[x,x,drop=FALSE],n1 =p$n1 ,n2 =p$n2 )$stat
    },p=param_i)

    matDist[i,]<-dist_i

    #if(shrink) setTxtProgressBar(pb, i)
    if(permute) pb$tick()
  }

  gdl<-sapply(index,function(x) { p<-length(x); p*(p+3)*0.5 })

  if(!shrink) shrink.param<-NULL

  return(list(matDist=matDist,gdl=gdl,param=ObsParam$shrink.info$lambda))
}

# Threshold for corrected p-value
#
# Identify the optimal cut-off
# @param mat statistic test distribution matrix
# @param alpha initial p-value threshold
# @param gdl ...
# @param permute ...
alphaCorrection<-function(mat,gdl,alpha=0.05,permute=FALSE){
  # convert pvalue

  if(permute){
    mat_pvalue<-
      apply(mat,2,function(col){
        sapply(col,function(x,cc) sum(cc>=x)/length(cc),cc=col )
      })
    type<-"minP"
  } else {
    mat_pvalue<-
      sapply(1:ncol(mat),function(i,m,g){
        1-stats::pchisq(m[,i],g[i])
      },m=mat,g=gdl)
    colnames(mat_pvalue)<-colnames(mat)
    type<-"maxT"
  }

  it<-0
  obs<-mat_pvalue[1,]
  names(obs)<-colnames(mat)
  ind_rej<-1:ncol(mat_pvalue)

  while( ncol(mat_pvalue)>0 & length(ind_rej)>0 ){
    it<-it+1
    dist_min<-apply(mat_pvalue,1,min)
    th_alpha<-stats::quantile(dist_min,alpha)
    ind_rej<-which(mat_pvalue[1,]<=th_alpha)
    if(length(ind_rej)>0) mat_pvalue<-mat_pvalue[,-ind_rej,drop=FALSE]
  }

  return(list(threshold=th_alpha,alpha=alpha,obs.pvalue=obs,iterations=it,type=type))
}

# Source Set for a single graph
#
# This function finds the source set of a given graph.
# @param ordering The output of the \code{\link[SourceSet]{ripAllRootsClique}} function.
# @param data a matrix of expression levels with column names for genes and row names for samples.
# @param classes a vector of 1 and 2 indicating the classes of samples. This vector must be matched with the rows of the data matrix, and can not contain more than two classes.
# @param seed integer value to get a reproducible random result.
# @param permute if set to TRUE, permutational p-values will be computed for the significance of the tests; otherwise, the asymptotic distribution will be used.
# @param theta copy from \code{\link[SourceSet]{sourceSet}}
# @param alpha the p-value threshold.
# @param shrink if set to TRUE, the algorithm will use the regularized estimate of the covariance matrices; otherwise, it will use the sample covariance matrices.
# @param shrink.function shrinkage function to be used if shrink=TRUE; default option is shrinkTEGS.
# @param shrink.param a list, with parameters to pass to shrink.function.
# @param print.label a string to use as name of the graph.
# @param return.permutations ....
# @details Refer to \code{\link[SourceSet]{sourceSet}}.
singleSourceSet<-function(ordering,data,classes,seed=NULL,
                          theta=1, permute=TRUE,
                          alpha=0.05,shrink=FALSE,shrink.function=shrinkTEGS, shrink.param=list(probs=0.05,type="min"),
                          print.label=NULL,return.permutations=FALSE){

  # max num of permutations
  max.perm<-10000

  #### check format data
  if(nrow(data)!=length(classes)) {
    stop("The row number must to be equal to length classes vector. Check data argument format.")
  }
  ### check all ordering variable in data (sistemare)
  if( ! all(ordering$elements$variables %in% colnames(data))  ){
    stop("All the variables in ordering argument must to be in data columns. Use graph::subGraph and run ripAllRootClique again.")
  }

  ## check allOrd
  if( !inherits(ordering,"allOrd") ){ #!(class(ordering)=="allOrd")
    stop("ordering argument is not an allOrd object. Use ripAllRootsClique function.")
  }
  ## check number of classes
  ncl<-table(classes)
  if( length(ncl)!=2 ){
    stop("the number of distinct classes must to be two")
  }

  # 2018-04-26 added drop
  dataR<-data[,colnames(data) %in%  ordering$elements$variables,drop=FALSE]

  ## Indices (exclude empty separator - "S0")
  n_comp<-length(ordering$indices$all)
  indcs<-lapply(ordering$indices$all[-n_comp],match,table=colnames(dataR))

  n<-nrow(dataR)
  p<-max(sapply(ordering$elements$cliques,length))

  ## Check if it is possible to use the asympotic test
  # Case: shrink=FALSE and nc<=p
  if( !shrink  & min(ncl)<= p ){
    msg<-paste0("The asymptotic mean-variance test can't be used. The number of samples in the smallest class (i.e. ",ncl,") must to be greater then number of variables in the maximal clique (i.e. ",p,"): shrink=TRUE option is used.")
    warning(msg)
    shrink<-TRUE
  }

  # Case shrink=TRUE & permute=FALSE
  if( shrink & !permute){
    warning("permute argument is not valid. The asymptotic distribution can't be used for the shrinkage estimate of the covariance matrices. permute=TRUE option is used.")
    permute<-TRUE
  }

  mtests<-nrow(ordering$elements$components)
  ## Set alpha correction option: T and type
  if(!permute) {
    #nperms.alpha<-nperms
    nperms.alpha<-round(min( max(500,(1/alpha)*(1+theta)) , max.perm))
    permute.alpha<-TRUE
    type.alpha<-"maxT"
  } else {
    #nperms.alpha<-nperms<-500
    nperms.alpha<- round(min(  max(1000,(mtests/alpha)*(1+theta)), max.perm  ))
    permute.alpha<-TRUE
    type.alpha<-"minP"
  }



  #### STEP 1: compute test statistic for each clique and each separator
  res<-testMeanVarianceDistribution(data = dataR,classes = classes, permute = permute.alpha ,indcs = indcs,nperms = nperms.alpha ,shrink = shrink,shrink.function = shrink.function, shrink.param = shrink.param, print.label = print.label )

  mat_test_element<-res$matDist
  gdl_test_element<-res$gdl

  # add S0 row (all test statistics equal to 0)
  mat_test_element<-cbind(mat_test_element,rep(0,nperms.alpha+1))
  colnames(mat_test_element)<-names(ordering$indices$all)
  gdl_test_element<-c(gdl_test_element,S0=0)

  #### STEP 2: combine clique/separator test statistics to obtain component statistics
  ## Test Statistic Matrix Component/Conditional: (nperm+1)x(ncomp)
  mat_test_component<-apply(ordering$elements$components,1,function(x,mat){
    ic<-as.numeric(x[3])
    is<-as.numeric(x[4])
    mat[,ic]-mat[,is]
  },mat=mat_test_element)
  colnames(mat_test_component)<-ordering$elements$components$component

  ## Gdl Vector Component/Conditional: (1)x(ncomp)
  gdl_test_component<-apply(ordering$elements$components,1,function(x,gdl){
    ic<-as.numeric(x[3])
    is<-as.numeric(x[4])
    gdl[ic]-gdl[is]
  },gdl=gdl_test_element)
  names(gdl_test_component)<-ordering$elements$components$component

  #### STEP 3: Indentify the corrected threshold for the p-value
  # if permute=FALSE it uses maxT
  # if permute=TRUE it uses minP
  ## Alpha correction
  corrected<-alphaCorrection(mat = mat_test_component,gdl = gdl_test_component, alpha=alpha,permute=permute)

  ## Resume All possible orderings + pvalue
  Decomposition<-ordering$indices$ordering
  for(i in 1:length(Decomposition)){
    ind<-match(Decomposition[[i]]$component,names(corrected$obs.pvalue))
    Decomposition[[i]]<-data.frame(Decomposition[[i]],pvalue=corrected$obs.pvalue[ind])
  }

  #### STEP 4: Union, take the variables contained in the cliques with significant p-value for each decomposition
  ## Union step
  union<-list()
  for(i in 1:length(Decomposition)){
    D<-Decomposition[[i]]

    ind_cl<-D$ind_clique[which(D$pvalue<=corrected$th)]
    ind_sep<-D$ind_separator[which(D$pvalue<=corrected$th)]

    u_i_cl<-unique(unlist(ordering$indices$all[ind_cl]))
    union<-append(union,list(u_i_cl))
    names(union)[i]<-names(Decomposition)[i]
  }

  #### STEP 5: Intersection, take the intersection of the sets identify in the previously step
  ## Intersection step
  sSet<-Reduce(intersect,union)
  if(length(sSet)==0) sSet<-NULL

  #### Extra: p-value graph
  dec<-ordering$indices$ordering[[1]]
  if(!permute){

    graph_gdl<-sum(gdl_test_component[dec$component])
    global_test_stat<-sum(mat_test_component[1,dec$component])
    PvalueGraph<-1-stats::pchisq(q = global_test_stat,df = graph_gdl )

  } else {
    # each ordering is a factorizationof the global distribution
    # i.e. the sum of the statistical test components is equal for each ordering
    global_test_stat_vec<-apply(mat_test_component[,dec$component,drop=FALSE],1,sum)
    PvalueGraph<-sum(global_test_stat_vec[-1]>=global_test_stat_vec[1])/length(global_test_stat_vec)
  }


  res<-list(primarySet=sSet,
            secondarySet=setdiff(unique(unlist(union)),sSet),
            orderingSet=union,
            Decompositions=Decomposition,
            Components=cbind(ordering$elements$components,pvalue=corrected$obs.pvalue),
            Elements=ordering$indices$all,
            Threshold=list(alpha=alpha,value=as.numeric(corrected$threshold),type=type.alpha,iterations=corrected$iterations,nperms=nperms.alpha),
            Graph=ordering$graph,
            PvalueGraph=PvalueGraph
  )

  if(return.permutations){
    res<-append(res,list(mat_test_element))
    names(res)[length(res)]<-"Permutations"
  }


  class(res)<-"sourceSet"
  return(res)
}


#' Source Set
#'
#' Identify the sets of variables that are potential sources of differential behavior,
#' (i.e., the primary genes) between two experimental conditions. The two experimental conditions are
#' associated to a set of graphs,
#' where each graph represents the topology of a biological pathway.
#' @param graphs a list of \code{graphNEL} objects  representing the pathways to be analyzed.
#' @param data  a matrix of expression levels with column names for genes and row names for samples; gene names must be unique.
#' @param classes a vector of length equal to the number of rows of \code{data}. It  indicates the class (condition) of each statistical unit. Only two classes, labeled as 1 and 2, are allowed;
#' @param seed integer value to get a reproducible random result. See \code{\link[base]{Random}}.
#' @param permute if \code{TRUE} permutation p-values are provided; if \code{FALSE}, asymptotic p-values are returned. NOTE: even if the argument permute is set to \code{FALSE} the function will permute the dataset; these permutations will be used to calculate the adjusted cut-off for the asymptotic p-values.
#' @param theta positive numeric value greater then 1, that defines the number of permutation. If \code{permute=TRUE}, (m/\code{alpha} x \code{theta}) permutations are used, where m is the number of unique conditional tests to be performed; otherwise, (1/\code{alpha} x \code{theta}) permutations are supplied.
#' @param alpha the p-value threshold. Denotes the level at which FWER is controlled for each input graph.
#' @param shrink if \code{TRUE}, regularized estimation of the covariance matrices is performed; otherwise, maximum likelihood estimations is used.
# @param shrink.function shrinkage function to be used if \code{shrink=TRUE}; default option is \code{\link[SourceSet]{shrinkTEGS}}.
# @param shrink.param a list, with parameters to pass to \code{shrink.function}.
#' @param return.permutations if \code{TRUE}, the function returns the matrix of test statistic values for the supplied (first row) and the permutated datasets.
#' @references
#' Sales, G. et al. (2017). graphite: GRAPH Interaction from pathway Topological Environment, r package version 1.22.0 edition.
#'
#' Westfall, P. and Young, S. (2017). Resampling-based multiple testing : examples and methods for p-value adjustment. Wiley.
#'
#' Djordjilovic, Vera and Chiogna, Monica (2022) Searching for a source of difference in graphical models. Journal of Multivariate Analysis 190, 104973
#'
#' Salviato, E. et al. (2019). \code{SourceSet}: a graphical model approach to identify primary genes in perturbed biological pathways. PLoS computational biology 15 (10), e1007357.
#' @details The \code{sourceSet} approach  models the data of the same pathway in two different
#' experimental conditions as realizations of two Gaussian graphical models sharing the same decomposable
#' graph G. Here, G = (V,E) is obtained from the pathway topology conversion, where V and E
#' represent genes and biochemical reactions, respectively.
#'
#' We give full freedom to the user in providing the underlying graph G, requiring only  a
#' specific input format (i.e., a \code{graphNEL} object). So, the user can provide a list of
#' manually curated pathways or use developed software to translate the bases of knowledge.
#' To date, the most complete software available for this task is \code{graphite} R package (Sales et al. 2017).
#'
#' The source set algorithm infers the set of primary genes (i.e., the source set) following - for each graph - five steps:
#' \itemize{
#'  \item{ decompose graph G in the set of the maximal cliques and the set of separators.}
#'  \item{ identify the cliques orderings, and the associated separators, that satisfy the running intersection property, using each cliques as root. See \code{\link[SourceSet]{ripAllRootsClique}}.}
#'  \item{ a) calculate marginal test statistics for the cliques and the separators, for both the original and the permutated datasets;
#'  b) compute the conditional test statistics for the unique components, calculated as the difference between clique and separator marginal test statistics;
#'  c) control the FWER, using the test statistics matrix of the previous point.}
#'  \item{ make the union of the sets of variables belonging to cliques that are associated to a significant test, within each decomposition. }
#'  \item{ derive the source set, defined as the intersection of the set of variables obtained in step 4 across decompositions.}
#' }
#'
#' Although the interpretation of the source set for a single graph is intuitive, the interpretation of the
#' collection of results associated to a set of pathways might be complex. For this reason,
#' we propose a guideline for the meta-analysis providing descriptive statistics and predefined plots. See,
#' \code{\link[SourceSet]{infoSource}}, \code{\link[SourceSet]{easyLookSource}}, \code{\link[SourceSet]{sourceSankeyDiagram}},  \code{\link[SourceSet]{sourceCytoscape}} and  \code{\link[SourceSet]{sourceUnionCytoscape}}.
#'
#'
#' @return The output of the function is an object of the \code{sourceSetList} class. It contains as many lists as the input graphs, and each of them provides the following variables:
#' \itemize{
#'  \item{\code{primarySet}: a character vector containing the names of the variables belonging to the estimated source set (primary dysregulation);}
#'  \item{\code{secondarySet}: a character vector containing the names of the variables belonging to the estimated secondary set (secondary dysregulation);}
#'  \item{\code{orderingSet}: a list of character vectors containing the names of the variables belonging to the estimated source set of each ordering; the union of these elements contains all genes affected by some form of perturbation; }
#'  \item{\code{Components}: a data frame that contains information about unique  tests, including their associated p-values;}
#'  \item{\code{Decompositions}: a list of data frames, one for each identified ordering. Each data frame is a subset of size k (i.e., number of cliques), of the \code{Components} elements}
#'  \item{\code{Elements}: cliques and separators of  the underlying decomposable graph. See \code{Graph}}
#'  \item{\code{Thresholds}: a list with information regarding the multiple testing correction:
#'     \itemize{
#'       \item{\code{alpha}: the input (nominal) significance level;}
#'       \item{\code{value}: the corrected threshold that ensures the control of FWER at level \code{alpha};}
#'       \item{\code{type}: the used procedure (minP or maxT);}
#'       \item{\code{iterations}: the number of iterations for the step-down procedure;}
#'       \item{\code{nperms}: the number of permutations.}
#'     }
#'  }
#'  \item{\code{Graph}: decomposable graph used in the analysis. It may differ from the input graph. In fact, if  the input graph is not  decomposable, the function will internally moralize and triangulate it.}
#'  }
#' @note If \code{permute} and/or \code{shrink} parameters violate the conditions required for the
#' existence of the full-rank maximum likelihood estimates, the algorithm reserves the possibility to change the user
#' settings through internal controls.
#'
#' Indeed, if the user wants to use the MLE of the covariance matrix (\code{shrink=FALSE}), all cliques -
#' in all pathways - must satisfy the \eqn{n > p_i} condition, where \eqn{n} is the number of samples for the
#' smaller class and \eqn{p_i} is the cardinality of the largest clique in the i-th pathway.
#' If even one clique does not satisfy this requirement, the regularized estimate must be used.
#' When a regularized estimate is employed (\code{shrink=TRUE}), the analytical null distribution
#' of the test statistics is no longer available, and we rely on permutation methods to obtain the
#' associated p-values.
#'
#' To address the multiple testing problem we use two versions of the method proposed by Westfall and Young (2017),
#' which uses permutations to obtain the joint distribution of the p-values.
#' More specifically, when the maximum likelihood estimates of the covariance matrices are used (\code{shrink=FALSE}), the asymptotic p-values and the maxT approach is adopted.
#' While, if the regularized estimates are calculated (\code{shrink=TRUE}), asymptotic distribution is no longer valid and the min P version and the
#' per-hypothesis permutation p-values to obtain the joint distribution of the p-values are needed.
#' The number of permutations depends on the method, the alpa level chosen, and the number of hypotheses. A minimum number of 500 and a maximum number of 10.000 permutations are allowed.
#'
#' @seealso \code{\link[graphite]{pathways}}, \code{\link[SourceSet]{infoSource}}, \code{\link[SourceSet]{easyLookSource}}, \code{\link[SourceSet]{sourceSankeyDiagram}},  \code{\link[SourceSet]{sourceCytoscape}} and  \code{\link[SourceSet]{sourceUnionCytoscape}}
#' @export
sourceSet<-function(graphs,data,classes,seed=NULL,
                    theta=1, permute=TRUE,
                    alpha=0.05,
                    shrink=FALSE,
                    return.permutations=FALSE){

  ## Check NA in data matrix
  if(sum(is.na(data))>0){
    stop("NAs not allowed in data matrix. Remove rows/columns that contain NAs or impute them (suggested function: impute.knn).")
  }

  ## Check columns names and nodes IDs
  g.id<-unique(unlist(sapply(graphs,function(x) graph::nodes(x))))
  n.mapped<- sum(colnames(data) %in% g.id)

  if( n.mapped==0 ) stop("Check data matrix columns names and graph nodes names, no match found")

  ## Default shrinkage Estimation of Covariance Matrix
  # Optimal choice of parameters for the source set analysis
  shrink.function<-shrinkTEGS
  shrink.param<-list(probs=0.05,type="min")

  start.time<-Sys.time()

  ## check graphs format
  if(!is.list(graphs)){
    stop("graphs argument must be a list of graphNEL object.")
    if( !inherits(graphs,"graphNEL") ){ # class(graphs)=="graphNEL"
     graphs<-list(graphs)
    }
  } else {
    N<-length(graphs)
    if(!all(sapply(graphs,class)=="graphNEL")){
      stop("all graphs list elements must to be a graphNEL object")
    }
  }

  ## add names to graph list elements
  if(is.null(names(graphs))){
    names(graphs)<-paste0("graph",1:N)
  }
  ind.na.name<-which(sapply(names(graphs),is.na))
  ind.empty.name<-which(sapply(names(graphs), function(x) x==""))
  ind.change<-sort(unique(c(ind.na.name,ind.empty.name)))
  if(length(ind.change)>0){
    names(graphs)[ind.change]<-paste0("graph",1:length(ind.change))
  }


  classes<-as.numeric(classes)
  ncl<-table(classes)
  ## check number of classes
  if( length(ncl)!=2 ){
    stop("the number of distinct classes must to be two")
  }

  ## check unique variable names (column)
  n.mult<-sum(table(colnames(data))>1)
  if(n.mult>0){
    stop("all variable names (column) must be unique")
  }

  ## STEP 0: Identify all possible orderings
  message("Identifing all possible orderings..")
  graphs<-subGraphData(data,graphs)

  ## check if (after sub graph) there are empty graphs
  ind.remove.graph<-which(sapply(sapply(graphs,function(x) graph::nodes(x)),length)==0)
  if(length(ind.remove.graph)>0){
    msg<-paste0("Removed ",length(ind.remove.graph)," graphs: empty graph or no match found between nodes and data matrix elements.")
    warning(msg)
    graphs<-graphs[-ind.remove.graph]
    if(length(graphs)==0){
      stop("all graphs are empty or there is no match between nodes and data matrix elements.")
    }
  }
  list_ordering<-lapply(graphs,ripAllRootsClique)


  ## check number of samples and number of parameters
  n<-nrow(data)
  pmax<-sapply(list_ordering,function(x){
    max(sapply(x$elements$cliques,length))
  })

  if(!shrink & min(ncl)<=max(pmax)){
    msg<-paste0("The asymptotic mean-variance test can't be used. The number of samples in the smallest class (i.e. ",min(ncl),") must to be greater then number of variables in the maximal clique (i.e. ",pmax,"): shrink=TRUE option is used.")
    warning(msg)
    shrink<-TRUE
  }


  ## Set number of permutations
  #mtest<-sapply(list_ordering,function(x) nrow(x$elements$components))

  res<-lapply(seq_along(list_ordering),function(i,lo){
    x<-lo[[i]]
    nx<-names(lo)[i]
    singleSourceSet(ordering = x,data = data,classes =classes ,seed = seed,theta =theta ,permute = permute,alpha = alpha,shrink = shrink,shrink.function = shrink.function,shrink.param = shrink.param,print.label = nx,return.permutations = return.permutations)
  },lo=list_ordering)
  names(res)<-names(list_ordering)

  end.time<-Sys.time()
  tot.time<-end.time-start.time
  msg<-paste("Total run time:",round(tot.time[[1]],2),attr(tot.time,"units"),sep=" ")
  message(msg)

  class(res)<-"sourceSetList"
  return(res)
}


# SubGraph
#
# Given the set of variables and a list of graphs create and returns subgraphs with the only supplied variables and any edges between them
# @param data A matrix or a dataframe with colnames for variables (genes) and rowname for samples
# @param graphs A list of graphs
# @seealso \code{\link[graph]{subGraph}}
# @export
subGraphData<-function(data,graphs){

  ## subdata and subgraph
  data_var<-colnames(data)
  reduced_graph<-lapply(graphs,function(x,dv) {
    gv<-graph::nodes(x)
    graph::subGraph(snodes = intersect(gv,dv),graph = x)
    },dv=data_var)

  return(reduced_graph)
}



#' Easy look results
#'
#' The function \code{easyLookSource}e allows to summarize the results obtained from the \code{\link[SourceSet]{sourceSet}} function through a heatmap, using \code{\link[ggplot2]{ggplot}} library.
#' @param sourceObj a \code{SourceSetObj} objects, i.e. the output of the \code{\link[SourceSet]{sourceSet}} function
#' @param name.graphs the graphs names to be visualized. Default value is \code{names(sourceObj)}
#' @param map.name.variable a list of customized labels to be associated with the names of the genes. Each list element must contain only one value (i.e. the new label), and the name of each element must be associated with the names of the genes given as input to the \code{\link[SourceSet]{sourceSet}} function (column names of \code{data} input argument). If a label is not mapped, the original name is used
#' @param label.variable title of the variable axis
#' @param label.graph title of the graph axis
#' @param maxnum.variable maximal number of variables to include in the plot. The variables are sorted internally and only the first \code{maxnum.variable} are be plotted.
#' @param maxnum.graph  maximal number of graphs to include in the plot. The graphs are sorted internally and only the first \code{maxnum.variable} are be plotted.
#' @param subname.variable number of characters of the variable names labels to show. The function cuts the name at the first \code{str.split} character that doesn't exceed \code{subname.variable}
#' @param subname.graph number of characters of the graph names labels to show. The function cuts the name at the first \code{str.split} character that doesn't exceed \code{subname.variable}
#' @param title overall title of the plot
#' @param subtitle subtitle of the plot
#' @param coord.equal if \code{TRUE}, forces the scale coordinate system to be equal for the y and x axis. See also \code{\link[ggplot2]{coord_fixed}}
#' @param coord.flip if \code{TRUE}, flips cartesian coordinates so that horizontal becomes vertical, and vertical, horizontal. Default option sets to x axis the variables, and to the y axis the graphs. See also \code{\link[ggplot2]{coord_flip}}
#' @param strsplit.graph character containing regular expression to use for cut graph labels to be shown. More details in \code{subname.variable} and  \code{subname.variable} argument descriptions
#' @param strsplit.variable character containing regular expression to use for cut varibale labels to be shown. More details in \code{subname.variable} and  \code{subname.variable} argument descriptions
#' @param col.primary cell color for the variables responsable of primary dysregulation
#' @param col.secondary cell color for the variables responsable of secondary dysregulation
#' @details
#' The plot is composed of a matrix whose rows represent pathways (i.e., graphs) and columns represent genes (i.e., variables).
#' Each cell i, j can take one of the following configurations:
#' \itemize{
#'   \item{ \code{2}: blue color, if the i-th gene is in the primary set of the j-th pathway}
#'   \item{ \code{1}: light blue color, if the i-th gene is in the secondary set of the j-th pathway}
#'   \item{ \code{0}: gray, if the i-th gene belongs to the j-th pathway }
#'   \item{ \code{NA}: white, if the i-th gene does not belong to the j-th pathway}
#' }
#' In the plot, the pathways are vertically ordered - top to bottom - according to the numbers of nodes in the source set. The genes are horizontally ordered (from left to right) based on the number of times they appear in a source set.
#' @seealso \code{\link[SourceSet]{sourceSet}}, \code{\link[SourceSet]{sourceSankeyDiagram}}
#' @return The function returns a \code{ggplot} object.
#' @examples
#' ## Load the SourceSetObj obtained from the source set analysis of ALL dataset
#'
#' # see vignette for more details
#' print(load(file=system.file("extdata","ALLsourceresult.RData",package = "SourceSet")))
#' class(results.all)
#' n.primary<-length(lapply(results.all,function(x) x$primarySet))
#'
#' # show only genes that appear in at least one of the source sets of the investigated pathways
#' easyLookSource(sourceObj=results.all, maxnum.variable = n.primary,
#'                label.variable = "Genes",label.graph = "Pathways")
#'
#' # flip coordinates
#' easyLookSource(sourceObj = results.all,maxnum.variable = n.primary,coord.flip = TRUE)
#' @export
easyLookSource<-function(sourceObj,name.graphs = names(sourceObj),
                         map.name.variable=NULL,
                         label.variable="Variable",label.graph="Graph",
                         subname.variable=10,subname.graph=20,
                         maxnum.variable=50, maxnum.graph=30,
                         title="Source Set for each Pathway",subtitle=NULL,
                         coord.equal=TRUE,coord.flip=FALSE,
                         strsplit.variable=" ",strsplit.graph=" ",
                         col.primary="#324E7B",col.secondary="#86A6DF"){

  if( !inherits(sourceObj,"sourceSetList") ){ #class(sourceObj)!="sourceSetList"
    stop("invalid class: sourceObj must be a sourceSetList object")
  }

  ### Check name.graphs
  ind.g<-which(name.graphs %in% names(sourceObj))
  ind.not<-setdiff(1:length(name.graphs),ind.g)
  if(length(ind.not)>0){
    msg<-paste0("One or more graph names are not in sourceObj (i.e., ",
                paste(name.graphs[ind.not],collapse = ", "),"): check the spelling")
  }
  sO<-sourceObj[name.graphs[ind.g]]

  if(length(sO)==0) stop(msg)
  else if(length(ind.not)>0) warning(msg)

  N<-length(sO)
  var_graph<-lapply(sO,function(x) unique(unlist(x$Elements)))
  var_data<-unique(unlist(var_graph))
  p<-length(var_data)

  source_matrix<-pvalue_matrix<-matrix(NA,ncol = p,nrow=N,dimnames = list(names(sO),var_data))

  for(i in 1:N){

    ## Fill source matrix:
    # NA= not in the graph; 1= in marginal set; 2= in source set; 0=  only in the graph
    si<-sO[[i]]$primarySet
    mi<-unique(unlist(sO[[i]]$orderingSet))

    source_matrix[i,var_graph[[i]]]<-0
    source_matrix[i,mi]<-1
    source_matrix[i,si]<-2

  }

  ## check if source and/or marginal sets are empty
  empty.source<-sum(source_matrix==2,na.rm=TRUE)==0
  empty.marginal<-sum(source_matrix==1,na.rm=TRUE)==0

  if( empty.marginal ){
    stop(paste0("No significant difference has been found between the two experimental conditions in all the ",N," analyzed graphs"))
  } else {
    if(empty.source) warning(paste0("The source sets are empty for all the ",N," analyzed graphs"))
  }


  info_variable<-data.frame(
    nsource=apply(source_matrix==2,2,sum,na.rm=TRUE),
    nmarginal=apply(source_matrix>0,2,sum,na.rm=TRUE),
    n=apply(!is.na(source_matrix),2,sum,na.rm=TRUE)
  )

  info_graph<-data.frame(
    nsource=apply(source_matrix==2,1,sum,na.rm=TRUE),
    nmarginal=apply(source_matrix>0,1,sum,na.rm=TRUE),
    n=apply(!is.na(source_matrix),1,sum,na.rm=TRUE)
  )

  ### Rank variables
  ind_variable<-order(info_variable$nsource,info_variable$nmarginal,info_variable$n,decreasing = TRUE)
  info_variable<-info_variable[ind_variable,]


  ## label variable
  fact_variable<-rownames(info_variable)


  ### Check map.name.variable (se non riesco a mappare un nome utilizzo di default 'none')
  if(!is.null(map.name.variable)){

    ind.name.var<-match(fact_variable,names(map.name.variable))
    ind.na<-which(is.na(ind.name.var))

    if(length(ind.na)>0){
      warn.msg<-paste0("One or more variable names are not mapped in map.name.variable: original name is used")
      warning(warn.msg)

      new.names<- paste0("id.",fact_variable[ind.na])
      names(new.names)<-fact_variable[ind.na]

      mpn<-append(map.name.variable,new.names)
      ind.name.var<-match(fact_variable,names(mpn))
    } else { mpn<-map.name.variable }

    mapped.name<-mpn[ind.name.var]
    mapped.name.vec<-sapply(mapped.name,paste,collapse=";")

  } else {
    mapped.name.vec<- fact_variable
    names(mapped.name.vec)<-fact_variable
  }


  lab_variable<-sapply(mapped.name.vec,function(x,nnn,sp){
    y<-strsplit(x,split = sp)[[1]]
    nch<-cumsum(sapply(y,nchar))
    nch<-nch+0:(length(nch)-1)

    ind<-which(nch<nnn | nch==nnn)
    if(length(ind)==0) st<-nch[[1]]
    else st<-as.numeric(nch[max(ind)])

    return(substr(x,start = 1,stop = st))
  },nnn=subname.variable,sp=strsplit.variable)


  ### Rank graphs
  ind_graph<-order(info_graph$nsource,info_graph$nmarginal,info_graph$n,decreasing = FALSE)
  info_graph<-info_graph[ind_graph,]

  fact_graph<-rownames(info_graph)
  lab_graph<-sapply(fact_graph,function(x,nnn,sp){
    y<-strsplit(x,split = sp)[[1]]
    nch<-cumsum(sapply(y,nchar))
    nch<-nch+0:(length(nch)-1)

    ind<-which(nch<nnn | nch==nnn)
    if(length(ind)==0) st<-nch[[1]]
    else st<-as.numeric(nch[max(ind)])

    return(substr(x,start = 1,stop = st))
  },nnn=subname.graph,sp=strsplit.graph)


  ## check label (duplicated label) ##
  name_dup<-names(which(table(lab_graph)>1))
  if(length(name_dup)>0){
    warning("label graph duplicated. To avoid this warning set an higher value for subname.graph")

    for(k in 1:length(name_dup)){
      nnk<-which(lab_graph %in% name_dup[k])
      lab_graph[nnk]<-paste(lab_graph[nnk]," (",1:length(nnk),")",sep="")
    }
  }

  name_dup<-names(which(table(lab_variable)>1))
  if(length(name_dup)>0){
    warning("label graph duplicated. To avoid this warning set an higher value for subname.graph")

    for(k in 1:length(name_dup)){
      nnk<-which(lab_variable %in% name_dup[k])
      lab_variable[nnk]<-paste(lab_variable[nnk],1:length(nnk))
    }
  }


  ## Reduce number of variables
  i_graph<-(nrow(source_matrix)-min(nrow(source_matrix),maxnum.graph)+1):(nrow(source_matrix))
  i_variable<-1:min(ncol(source_matrix),maxnum.variable)


  source_matrix_cut<-source_matrix[rownames(info_graph)[i_graph],
                                   rownames(info_variable)[i_variable],drop=FALSE]

  lab_variable<-lab_variable[names(lab_variable) %in% colnames(source_matrix_cut)]
  lab_graph<-lab_graph[names(lab_graph) %in% rownames(source_matrix_cut)]

  dataframe_source<-reshape2::melt(source_matrix_cut)
  colnames(dataframe_source)<-c("Graph","Variable","Source")


  ## coord_equal and coord_flip are exclusive
  if(!coord.flip){
    dataframe_source$Variable<-factor(x = dataframe_source$Variable,levels =names(lab_variable),labels = lab_variable )
    dataframe_source$Graph<-factor(x = dataframe_source$Graph,levels =names(lab_graph), labels = lab_graph )

    source_plot<-ggplot2::ggplot(data = dataframe_source,ggplot2::aes_(x=~as.factor(Variable),y=~as.factor(Graph)))+
      ggplot2::geom_tile(ggplot2::aes_(fill=~as.factor(Source)),colour="white")+
      ggplot2::scale_fill_manual(name="Source set",limits=c("0","1","2"),values=c("gray",col.secondary,col.primary),na.value = "white")+
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "gray"),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,vjust = 0)  )+
      ggplot2::xlab(label.variable)+ggplot2::ylab(label.graph)+ggplot2::labs(label=title,subtitle = subtitle)
  } else{
    dataframe_source$Variable<-factor(x = dataframe_source$Variable,levels =rev(names(lab_variable)),labels = rev(lab_variable) )
    dataframe_source$Graph<-factor(x = dataframe_source$Graph,levels =rev(names(lab_graph)), labels = rev(lab_graph) )

    source_plot<-ggplot2::ggplot(data = dataframe_source,ggplot2::aes_(y=~as.factor(Variable),x=~as.factor(Graph)))+
      ggplot2::geom_tile(ggplot2::aes_(fill=~as.factor(Source)),colour="white")+
      ggplot2::scale_fill_manual(name="Source set",limits=c("0","1","2"),values=c("gray",col.secondary,col.primary),na.value = "white")+
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "gray"),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,vjust = 0)  )+
      ggplot2::xlab(label.graph)+ggplot2::ylab(label.variable)+ggplot2::labs(label=title,subtitle = subtitle)
  }


  if(coord.equal) source_plot<-source_plot+ggplot2::coord_equal()



  return(source_plot)
}


#' Create a D3 JavaScript Sankey diagram
#'
#' The function \code{sourceSankeyDiagram} allows to summarize the results obtained from the \code{\link[SourceSet]{sourceSet}} function through a Sankey diagram, highlighting the relationships among nodes, graphs, and source sets.
#' @param sourceObj a \code{SourceSetObj} objects, i.e. the output of the \code{\link[SourceSet]{sourceSet}} function
#' @param name.graphs the names of the graphs to be visualized. Default value is \code{names(sourceObj)}
#' @param map.name.variable a list of customized labels to be associated with the names of the genes. Each list element must contain only one value (i.e. the new label), and the name of each element must be associated with the names of the genes given as input to the \code{\link[SourceSet]{sourceSet}} function (column names of \code{data} input argument). If a label is not mapped, the original name is used
#' @param cutoff the maximum number of variables to include in the sankey graph. The final number of visualized variables could be higher than the \code{cutoff} number
#' @param cut.extra.module if set to \code{TRUE}, modules consisting only of variables excluded by the cutoff are not displayed
#' @param height numeric height (in pixels) for the network graph's frame area
#' @param width numeric width (in pixels) for the network graph's frame area
#' @seealso \code{\link[networkD3]{sankeyNetwork}}, \code{\link[SourceSet]{sourceSet}}, \code{\link[SourceSet]{easyLookSource}}
#' @details
#' The layout is organized on three levels:
#' \itemize{
#'    \item{ the first level (left) shows nodes that appear in at least one source sets of the analyzed graphs; }
#'    \item{ the second level (center) is made up of modules. A module is defined as a set of nodes belonging to a connected subgraph of one pathway, that is also contained in associated source set. A pathway can have multiple modules, and, at the same time, one module can be contained in multiple pathways;}
#'    \item{ the third level (right) shows of pathways. }
#'  }
#'  The three levels are to be read from left to right. A link between left element a and right element b must be interpret as "a is contained in b".
#'
#'  The implementation of the \code{sourceSankeyDiagram} function takes advantage of the D3 library (JavaScript), making the plot interactive.
#'  In fact, it is possible to vertically shift the displayed elements, and to view some useful information by positioning the cursor over items and links.
#' @references
#' Allaire, J.J., Gandrud, G., Russell, K., and Yetman, C.J. (2017). networkD3: D3 JavaScript Network Graphs from R, r package version 0.4 edition.
#'
#' Bostock, M., Ogievetsky, V., and Heer, J. (2011). D3 data-driven documents. IEEE Transactions on Visualization and Computer Graphics, 17(12):2301–2309.
#' @examples
#' ## Load the SourceSetObj obtained from the source set analysis of ALL dataset
#'
#' # see vignette for more details
#' print(load(file=system.file("extdata","ALLsourceresult.RData",package = "SourceSet")))
#' class(results.all)
#'
#' sourceSankeyDiagram(sourceObj = results.all ,cut.extra.module = FALSE )
#'
#' # shows the variable that appears most often in the source sets
#' sourceSankeyDiagram(sourceObj = results.all, cutoff = 1 ,cut.extra.module = FALSE )
#' # cut modules in which the variable is not contained
#' sourceSankeyDiagram(sourceObj = results.all, cutoff = 1 ,cut.extra.module = TRUE )
#' @return The function returns an interactive \code{sankeyNetwork} object.
#' @export
sourceSankeyDiagram<-function(sourceObj,name.graphs=names(sourceObj),
                              map.name.variable=NULL,
                              cutoff=50,
                              cut.extra.module=TRUE,
                              height=NULL,width=NULL){

  if( !inherits(sourceObj,"sourceSetList") ){ #class(sourceObj)!="sourceSetList"
    stop("invalid class: sourceObj must be a sourceSetList object")
  }

  ## hight constant
  cH<-50

  ### Check name.graphs
  ind.g<-which(name.graphs %in% names(sourceObj))
  ind.not<-setdiff(1:length(name.graphs),ind.g)
  if(length(ind.not)>0){
    msg<-paste0("One or more graph names are not in sourceObj (i.e., ",
                paste(name.graphs[ind.not],collapse = ", "),"): check the spelling")
  }

  sO<-sourceObj[name.graphs[ind.g]]

  if(length(sO)==0) stop(msg)
  else if(length(ind.not)>0) warning(msg)

  #### UPDATE: suggestions
  # cancellare layer module tra le opzioni
  # fornire una lista di subset (variabili o grafi)

  # source set
  N<-length(sO)
  S<-lapply(sO,function(x) x$primarySet)
  S<-S[!sapply(S,is.null)]


  ## check if source sets are empty
  empty.source<-length(S)==0
  if( empty.source ){
    stop(paste0("The source sets are empty for all the ",N," analyzed graphs"))
  }

  # graph (with a non null source set)
  G<-names(S)
  # variable (with a non null source set)
  V<-unique(unlist(S))


  ## cut the list, to obtain a interpretable graph
  if(length(V)>cutoff){
    cut<-sort(table(unlist(S)),decreasing = TRUE)
    nco<-as.numeric(cut[min(length(V),cutoff)])
    ind_co<-max(cutoff,max(which(cut==nco)) )
    var_co<-names(cut[1:ind_co])
  } else { var_co<- V }

  N<-length(G)
  modules<-list()

  for(i in 1:N){


    if( !is.null(S[[i]]) ){
      subg<-graph::subGraph(S[[i]],sO[[names(S)[i]]]$Graph)
      cl<-igraph::clusters(igraph::igraph.from.graphNEL(subg))

      nm<-cl$no
      mi<-lapply(1:nm,function(i,el){
        sort(names(which(el==i)))
      },el=cl$membership)

    } else mi<-list(NULL)
    modules<-append(modules,list(mi))
    names(modules)[i]<-G[i]

  }

  ### STEP 2: Identify unique module
  # unique funziona perchè li ho ordinato per nome gli elementi in lapply
  uM<-unlist(modules,recursive = FALSE)
  uM<-unique(uM[!sapply(uM,is.null)])
  names(uM)<-paste0("module",1:length(uM))


  # Convert variable names
  #map.name.variable=NULL,
  ### Check map.name.variable (se non riesco a mappare un nome utilizzo di default 'none')
  if(!is.null(map.name.variable)){

    map.uM<-unique(unlist(uM))
    ind.name.var<-match(map.uM,names(map.name.variable))
    ind.na<-which(is.na(ind.name.var))

    if(length(ind.na)>0){
      warn.msg<-paste0("One or more variable names are not mapped in map.name.variable: original name is used")
      warning(warn.msg)

      new.names<- paste0("id.",map.uM[ind.na])
      names(new.names)<-map.uM[ind.na]

      mpn<-append(map.name.variable,new.names)
      ind.name.var<-match(map.uM,names(mpn))
    } else { mpn<-map.name.variable }

    mapped.name<-mpn[ind.name.var]
    uM<-lapply(uM,function(x) unlist(mpn[x]))
    V<-unlist(mpn[V])
    modules<-lapply(modules,function(mod.path){
        lapply(mod.path,function(x) unlist(mpn[x]) )
      })
    var_co<-unlist(mpn[var_co])

  }


  ### STEP 3: build dataset
  # identify links between variables (V) and modules (uM)
  link_var_mod<-lapply(seq_along(uM),function(i,vv,mm){
    ind<-which(vv %in% mm[[i]])
    cbind(vv[ind],names(mm)[i])
  },vv=V,mm=uM)

  # identify links between modules (uM) and graph (G)
  # devo usare l'oggetto module
  link_mod_graph<-lapply(seq_along(modules),function(i,gg,mm){

    name<-names(gg)[i]
    x<-sapply(gg[[i]],paste,collapse=";")
    ind<-match(x,mm)
    cbind(names(mm)[ind],rep(name,length(ind)))

  },gg=modules,mm=sapply(uM,paste,collapse=";"))

  sankey_data<-data.frame(Reduce("rbind",append(link_var_mod,link_mod_graph)),stringsAsFactors = FALSE,row.names = NULL)
  colnames(sankey_data)<-c("from","to")

  ## CUT: mantain only the element selected by the cutoff
  link_var_mod_cut<-lapply(link_var_mod,function(x,lco){
    ind<-which(x[,1] %in% lco)
    if(length(ind)>0) return( x[ind,,drop=FALSE])
    else return(NULL)
  },lco=var_co )
  link_var_mod_cut<-link_var_mod_cut[!sapply(link_var_mod_cut,is.null)]

  if(cut.extra.module){

    mod_co<-unique(unlist(lapply(link_var_mod_cut,function(x){ x[,2] })))
    link_mod_graph_cut<-lapply(link_mod_graph,function(x,lco){
      ind<-which(x[,1] %in% lco)
      if(length(ind)>0) return( x[ind,,drop=FALSE])
    },lco=mod_co)
    link_mod_graph_cut<-link_mod_graph_cut[!sapply(link_mod_graph_cut,is.null)]

  } else  { link_mod_graph_cut<-link_mod_graph }

  sankey_data_cut<-data.frame(Reduce("rbind",append(link_var_mod_cut,link_mod_graph_cut)),stringsAsFactors = FALSE,row.names = NULL)
  colnames(sankey_data_cut)<-c("from","to")

  ### SANKEY: network D3
  # peso di ciascun nodo= pesonodo
  # peso di ciascun modulo = length(modulo) * pesonodo
  # pedi di ciascun grafo = sum(length(moduli)) * pesonodo
  SNode<-data.frame(name=unique(c(sankey_data_cut$from,sankey_data_cut$to)))
  nw<-rep(1,length(V)); names(nw)<-V


  nm<-sapply(uM,length)
  tt<-table(sankey_data_cut$from)
  tt<-tt[grep("module",names(tt))]
  nm<-nm[names(tt)]
  #tt<-tt[match(names(uM),names(tt))]

  #weight<-c(nw,sapply(uM,length),sapply(S,length))
  weight<-c(nw,nm/tt,sapply(S,length))

  ## height setting
  if(is.null(height)){
    height<-max(length(uM),length(G),length(var_co))*cH
  }



  SLink<-data.frame(
    source=match(sankey_data_cut$from,SNode$name)-1,
    target=match(sankey_data_cut$to,SNode$name)-1,
    value=as.numeric(weight[match(sankey_data_cut$from,names(weight))])
    #value= as.numeric(weight[match(SNode$name,names(weight))])
  )

  networkD3::sankeyNetwork(Links =SLink ,Nodes = SNode, Source = "source",
                Target = "target", Value = "value", NodeID = "name",
                units = "", fontSize = 10, nodeWidth = 10,
                height = height, width = width)

  #return(
  #  list(sankey=sankey_data,modules=uM)
  #)

}


# Moralize a graph
#
# Moralize ....
# @param graph ...
moralize<-function(graph){

  if( !inherits(graph,"graphNEL") ){ #class(graph)!="graphNEL"
    stop("graph argument is not a graphNEL object")
  }

  if(gRbase::is.UG(graph)){
    #return(graph)
    return(gRbase::M2graphNEL(gRbase::as.adjMAT(graph)))
  }

  if(!gRbase::is.DAG(graph)){
    m <- gRbase::graphNEL2M(graph)
    diag(m)<-0

    n<-nrow(m)
    m.dag<-m
    m.ug<-matrix(0,ncol=n,nrow=n,dimnames = list(rownames(m),colnames(m)))

    for(i in 1:(n-1)){
      for(j in (i+1):n){
        if( m[i,j]!=0 &  m[i,j]==m[j,i]){
          m.ug[i,j]<-m.ug[j,i]<-1
          m.dag[i,j]<-m.dag[j,i]<-0
        }
      }
    }


    m.moralize<-gRbase::moralizeMAT(m.dag)+m.ug
    g.moralize<-methods::as(m.moralize,"graphNEL")

  } else {
    g.moralize<-gRbase::moralize(graph)
  }

  return(g.moralize)

}


# Plot the source graph
#
# This function allow to easily highligth the source and the marginal set of your graph
# @param sourceObj ...
# @param name the name of the graph to plot in the sourceObj list. It must be contained in names(sourceObj).
# @param maxnum.variable the maximun number of nodes to plot.
# @param type type of graph. It can be \code{igraph} or \code{graphNEL}.
# @param dag.graph optional graph. In the sourceObj list is contained the triangulated graph. If you want to use the original graph use this argument. The argument must to have the same nodes of the triangualted graph of the \code{sourceObj}.
# @param graph.size graph size
sourceGraph<-function(sourceObj,name=NULL,
                      maxnum.variable=50,
                      graph.size=7,
                      type="igraph",
                      dag.graph=NULL){

  if(  !inherits(sourceObj,"sourceSetList")  ){ #class(sourceObj)!="sourceSetList"
    stop("invalid class: sourceObj must be a sourceSetList object")
  }

  ## table colors ##
  ind<- which(names(sourceObj)==name)
  if(length(ind)!=1){
    stop("invalid argument: name must to be contained in names(sourceObj)")
  }


  if(!is.null(dag.graph)){
    if( !inherits(dag.graph,"graphNEL")  ){ #class(dag.graph)=="graphNEL"
      graph<-dag.graph
    } else {
      stop("dag.graph must to be a graphNEL object")
    }
  } else {
    graph<-sourceObj[[ind]]$Graph
  }

  sO<-sourceObj[[ind]]
  marginal<-table(unlist(sO$orderingSet))

  ## Color palette
  CRPal <- grDevices::colorRampPalette(c("#86A6DF","#5068A9"))
  COL<-c(CRPal(max(marginal)-1),"#324E7B")

  N<-graph::nodes(graph)
  #### node attribute ###
  nA<-list()
  # node fillcolor: number of times that the node appears in the marginal set
  nA$fillcolor<-stats::setNames(rep("white",length(N)),N)
  nA$fillcolor[names(marginal)]<-COL[marginal]

  # node font color:
  nA$fontcolor<-stats::setNames(rep("#9A9B94",length(N)),N)
  nA$fontcolor[names(marginal)]<-"white"

  # node size
  nA$height<-stats::setNames(rep(0.5,length(N)),N)
  nA$height[names(which(marginal==max(marginal)))]<-0.8

  #### general attribute ###
  A<-list()
  A$edge$color<-"#9A9B94"
  A$node$color<-"#324E7B"
  A$graph$size<-c(graph.size,graph.size)

  Rgraphviz::plot(graph,nodeAttrs=nA,attrs=A)

}


#' Get summary statistics on graphs and variables
#'
#' The \code{infoSource} function provides a summary of the results by focusing on either variables or graphs.
#' @param sourceObj a \code{SourceSetObj} object, i.e. the output of the \code{\link[SourceSet]{sourceSet}} function
#' @param map.name.variable a list of customized labels to be associated with the names of the genes. Each list element must contain only one value (i.e. the new label), and the name of each element must be associated with the names of the genes given as input to the \code{\link[SourceSet]{sourceSet}} function (column names of \code{data} input argument). If a label is not mapped, the original name is used
#' @param method correction method for p-values calculated on graphs. The adjustment methods allowed are: \code{fdr} (default), \code{holm}, \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH}, \code{BY} or \code{none}. For more details refer to \code{\link[stats]{p.adjust}}.
#' @return
#' The function guides the user in identifying interesting variables returning two objects:
#' \itemize{
#'   \item{\code{graph}: a dataframe that summirizes the results of the individual input graphs, composed as follows:
#'     \itemize{
#'       \item{\code{n.primary}: number of genes belonging to the source set;}
#'       \item{\code{n.secondary}: number of genes belonging to the secondary set;}
#'       \item{\code{n.graph}: number of genes within the graph;}
#'       \item{\code{n.cluster}: number of connected components of the graph;}
#'       \item{\code{primary.impact}: relative size of the estimated source set. This index quantifies the proportion of the graph impacted by primary dysregulation;}
#'       \item{\code{total.impact}: relative size of the set of genes impacted by dysregulation. This index quantifies the proportion of the graph impacted by either primary or secondary dysregulation;}
#'       \item{\code{adj.pvalue}: multiplicity adjusted p-value for the hypothesis of equality of the two distributions associated to the given graph}
#'     }}
#'   \item{\code{variable}: a dataframe that summarized the results of the individual variables, composed as follows:
#'     \itemize{
#'       \item{\code{n.primary}: number of input graphs in which the gene appears in the associated source set;}
#'       \item{\code{n.secondary}: number of input graphs in which the gene appears in the associated secondary set;}
#'       \item{\code{n.graph}:  number of pathways in which the gene is annotated;}
#'       \item{\code{specificity}: percentage of input graphs containing the given genes with respect to the total number of input graphs;}
#'       \item{\code{primary.impact}: percentage of input graphs, such that the given gene belongs to their estimated source set, with respect to the total number of input graphs in which the gene appears;}
#'       \item{\code{total.impact}: percentage of input graphs, such that the given gene is affected by some form of dysregulation in the considered graph, with respect to the total number of input graphs in which the gene appears;}
#'       \item{\code{relevance}: percentage of the input graphs such that the given variable belongs to their estimated source set, with respect to the total number of input graphs. It is a general measure of the importance of the gene based on the chosen pathways;}
#'       \item{\code{score}: a number ranging from \code{0} (low significance) to \code{+Inf} (maximal significance), computed as the combination of the p-values of all components (of all the input graphs) containing the given variable}
#'     }
#'   }
#' }
#' @note
#' Ideally, variables of the primary dysregulation will be elements of the source set in all input graphs that contain them and will thus have high values of \code{source.impact} and \code{score}.
#' However, if a given variable appears in a single graph, and belongs to its source set, these indices can be deceptive.
#'
#' For this reason, \code{relevance} serves to identify variables that apart from being good candidates for primary genes, also appear frequently in the input graphs.
#' Which index is to be preferred depends on the objective of the analysis: in case of exploratory analysis, we suggest to rely on \code{relevance}.
#' @examples
#' ## Load the SourceSetObj obtained from the source set analysis of ALL dataset
#'
#' # see vignette for more details
#' print(load(file=system.file("extdata","ALLsourceresult.RData",package = "SourceSet")))
#' class(results.all)
#'
#' info.all<-infoSource(sourceObj = results.all)
#' ## results of individual input graphs
#' info.all$graph
#'
#' ## results of individual variables
#' # ..that appear in more than one graph and with relevance>0
#' info.all.genes<-info.all$variable[info.all$variable$n.graph>1 & info.all$variable$relevance>0,]
#' # ..ordered by score
#' ind.ord<-order(info.all.genes$relevance,decreasing = TRUE)
#' info.all.genes[ind.ord,]
#' @export
infoSource<-function(sourceObj,map.name.variable=NULL,method="fdr"){

  if( !inherits(sourceObj,"sourceSetList") ){ #class(sourceObj)!="sourceSetList"
    stop("invalid class: sourceObj must be a sourceSetList object")
  }

  N<-length(sourceObj)
  var_graph<-lapply(sourceObj,function(x) unique(unlist(x$Elements)))
  var_data<-unique(unlist(var_graph))
  p<-length(var_data)

  source_matrix<-pvalue_matrix<-matrix(NA,ncol = p,nrow=N,dimnames = list(names(sourceObj),var_data))

  for(i in 1:N){

    ## Fill source matrix:
    # NA= not in the graph; 1= in marginal set; 2= in source set; 0=  only in the graph
    si<-sourceObj[[i]]$primarySet
    #mi<-unique(unlist(sourceObj[[i]]$orderingSet))
    mi<-sourceObj[[i]]$secondarySet

    source_matrix[i,var_graph[[i]]]<-0
    source_matrix[i,mi]<-1
    source_matrix[i,si]<-2

  }

  nclust<-sapply(sourceObj,function(x){
    cl<-igraph::clusters(igraph::igraph.from.graphNEL(x$Graph))
    length(cl$csize)
  })



  ## INFO: variables & graph ##
  info_variable<-data.frame(
    n.primary=apply(source_matrix==2,2,sum,na.rm=TRUE),
    n.secondary=apply(source_matrix==1,2,sum,na.rm=TRUE),
    n.graph=apply(!is.na(source_matrix),2,sum,na.rm=TRUE)
  )


  info_graph<-data.frame(
    n.primary=apply(source_matrix==2,1,sum,na.rm=TRUE),
    n.secondary=apply(source_matrix==1,1,sum,na.rm=TRUE),
    n.graph=apply(!is.na(source_matrix),1,sum,na.rm=TRUE),
    n.cluster=nclust
  )

  ## pvalue/score ##
  var_pval<-scoreNode(sourceObj)

  #gr_pval<-pvalueGraph(sourceObj = sourceObj,method=method)
  gr_pval<-stats::p.adjust(sapply(sourceObj,function(x) x$PvalueGraph),method = method)

  ### Check map.name.variable (se non riesco a mappare un nome utilizzo di default 'none')
  if(!is.null(map.name.variable)){

    ind.name.var<-match(rownames(info_variable),names(map.name.variable))
    ind.na<-which(is.na(ind.name.var))

    if(length(ind.na)>0){
      warn.msg<-paste0("One or more variable names are not mapped in map.name.variable: original name is used")
      warning(warn.msg)

      new.names<- paste0("id.",rownames(info_variable)[ind.na])
      names(new.names)<-rownames(info_variable)[ind.na]

      mpn<-append(map.name.variable,new.names)
      ind.name.var<-match(rownames(info_variable),names(mpn))
    } else { mpn<-map.name.variable }

    mapped.name<-mpn[ind.name.var]
    mapped.name.vec<-sapply(mapped.name,paste,collapse=";")

  } else {
    mapped.name.vec<- rownames(info_variable)
    names(mapped.name.vec)<- rownames(info_variable)
  }

  # add impact index
  info_variable<-data.frame(
    mapped.name=mapped.name.vec,
    info_variable,
    specificity=round(info_variable$n.graph/length(sourceObj),5),
    primary.impact=round(info_variable$n.primary/info_variable$n.graph,5),
    total.impact=round((info_variable$n.primary+info_variable$n.secondary)/info_variable$n.graph,5),
    score=round(var_pval$score$score[match(rownames(info_variable),rownames(var_pval$score))],5)
  )
  info_variable<-data.frame(
    info_variable,
    relevance=round(info_variable$primary.impact*info_variable$specificity,5)
  )



  # add impact index
  info_graph<-data.frame(
    info_graph,
    primary.impact=round(info_graph$n.primary/info_graph$n.graph,5),
    total.impact=round((info_graph$n.primary+info_graph$n.secondary)/info_graph$n.graph,5),
    adj.pvalue=gr_pval[rownames(info_graph)]
  )






  ## Rank variables
  #ind_variable<-order(info_variable$specificity,info_variable$source.impact,info_variable$score,decreasing = TRUE)
  ind_variable<-order(info_variable$relevance,info_variable$score,info_variable$n.graph,decreasing = TRUE)
  info_variable<-info_variable[ind_variable,]



  list(variable=info_variable,graph=info_graph)
}


# Score nodes
#
# ...
# @param sourceObj ..
scoreNode<-function(sourceObj){

  if( !inherits(sourceObj,"sourceSetList") ){ #class(sourceObj)!="sourceSetList"
    stop("invalid class: sourceObj must be a sourceSetList object")
  }

  transformation<-function(x,alpha,th){
    if(x<=th){
      y<-scales::rescale(x,from=c(0,th),to=c(0,alpha))
    } else {
      y<-scales::rescale(x,from=c(th,1),to=c(alpha,1))
    }
    y
  }


  N<-length(sourceObj)
  var_graph<-lapply(sourceObj,function(x) unique(unlist(x$Elements)))
  var_data<-unique(unlist(var_graph))
  p<-length(var_data)

  pvalue_matrix<-matrix(NA,ncol = p,nrow=N,dimnames = list(names(sourceObj),var_data))

  for(i in 1:N){

    alpha_i<-sourceObj[[i]]$Threshold$alpha
    theta_i<-sourceObj[[i]]$Threshold$value

    ## Fill pvalue matrix
    # for each decomposition assign the pvalue of the component at each variable (it appears only once)
    di<-sourceObj[[i]]$Decompositions
    compi<-sourceObj[[i]]$Components
    K<-length(di)

    ## create the set for each component (clique.el-sep.el) and assign the same pvalue
    pi<-plyr::alply(compi,1,function(x,el){
      cl<-as.numeric(x[[3]]); sep<-as.numeric(x[[4]])
      set<-setdiff(el[[cl]],el[[sep]])
      pval<-rep(as.numeric(x[[6]]),length(set))
      pval<-as.vector(pval)
      names(pval)<-set
      pval
    },el=sourceObj[[i]]$Elements)
    names(pi)<-compi$component

    ## assign at each decomposition and each variable the pvalue
    pval_decomp<-matrix(NA,ncol=length(var_graph[[i]]),nrow = K,dimnames = list(paste0("dec",1:K),var_graph[[i]]))
    for(k in 1:K){
      pvalk<-Reduce("c",pi[di[[k]]$component])
      pval_decomp[k,names(pvalk)]<-pvalk
    }

    ## assign at each variable the max pvalue obtained in the K decomposition
    pval_graph<-apply(pval_decomp,2,max)

    ## transform p-value on a comparable scale
    pval_graph_star<-sapply(pval_graph,transformation,alpha=alpha_i,th=theta_i)


    pvalue_matrix[i,names(pval_graph_star)]<-pval_graph_star

  }

  resume<-data.frame(
    score= -log(apply(pvalue_matrix,2,mean,na.rm=TRUE)),
    frequency=apply(pvalue_matrix,2,function(col) sum(!is.na(col)))
  )

  th<- -log(mean(sapply(sourceObj,function(x) x$Threshold$value)))
  #resume<-resume[resume$score>=th,]
  #resume<-resume[order(-resume$frequency,resume$score),]

  list(score=resume,meanth=th)
}


#' Visualize in Cytoscape a collection of graphs analyzed with the source set algorithm
#'
#' The function, thanks to the connection with the Cytoscape software, allows the user to create a collection of graphs to be visualized in a unique session, while documenting interesting findings.
#' @param sourceObj a \code{SourceSetObj} objects, i.e. the output of the \code{\link[SourceSet]{sourceSet}} function.
#' @param name.graphs the names of the graphs to be visualized. Default value is \code{names(sourceObj)}. NOTE: even if a subset of graphs are selected in \code{name.graphs}, the returned statistics are always calculated on the entire collection in the \code{sourceObj} argument.
#' @param map.name.variable a list of customized labels to be associated with the names of the genes. Each list element must contain only one value (i.e. the new label), and the name of each element must be associated with the names of the genes given as input to the \code{\link[SourceSet]{sourceSet}} function (column names of \code{data} input argument). If a label is not mapped, the original name is used.
#' @param collection.name name of the collection of graphs displayed in Cytoscape.
#' @param method correction method for p-values calculated on graphs. The adjustment methods allowed are: \code{fdr} (default), \code{holm}, \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH}, \code{BY} or \code{none}. For more details refer to \code{\link[stats]{p.adjust}}.
#' @details The visual node attributes size and fill color are defined in a dynamic manner through a visual mapping based on the indices provided by the \code{\link[SourceSet]{infoSource}} function (automatically uploaded in the bottom panel - right side).
#'
#' A discrete mapper between \code{source} attribute and size is applied:
#' \itemize{
#'   \item{big size: the variable belongs to the primary set (\code{source=2});}
#'   \item{medium size: the variable belongs to the secondary set (\code{source=1});}
#'   \item{small size: otherwise (\code{source=0}).}
#' }
#' On the other hand, a color gradient mapper between fill node color and \code{relevance} is adopted: higher values are highlighted with darker blue color.
#'
#' The default style can be changed manually either within Cytoscape (for further information see \href{http://manual.cytoscape.org/en/stable/Styles.html}{manual}) or within an R package \code{r2cytoscape} through network SUID returned by the \code{sourceCytoscape function} (for further details see \href{https://github.com/cytoscape/r2cytoscape}{manual}).
#'
#' It is also possible to call the sourceCytoscape function multiple times, with all the graphs being visualized in a unique session within a collection specified by collection.name.
#' @note The function use the \code{r2cytoscape} package to connect to Cytoscape from R using CyREST. \code{r2cytoscape} can be downloaded from:
#' \itemize{
#'    \item{Bioconductor: \code{biocLite("r2cytoscape")};}
#'    \item{GitHub: \code{install_github("cytoscape/r2cytoscape")}.}
#' }
#'
#' To enable the display function to work properly, three simple steps are required:
#' \itemize{
#'    \item{ Download \href{https://cytoscape.org/download.php}{Cytoscape} (version 3.3 or later);}
#'    \item{ Complete installation wizard;}
#'    \item{ Launch Cytoscape (before calling the functions).}
#' }
#' @seealso \code{\link[SourceSet]{sourceSet}}, \code{\link[SourceSet]{infoSource}}. \code{\link[SourceSet]{sourceUnionCytoscape}}, \code{r2cytoscape}
#' @return The function returns an interactive session in Cytoscape.?
#' @examples
#' ## Load the SourceSetObj obtained from the source set analysis of ALL dataset
#'
#' # see vignette for more details
#' print(load(file=system.file("extdata","ALLsourceresult.RData",package = "SourceSet")))
#' class(results.all)
#'
#' ## NB: Remember to launch cytoscape before running the following commands
#' # Create two collections of pathways to visualize the results
#' graph.signaling<-names(results.all)[grep("signaling",names(results.all))]
#' graph.other<-setdiff(names(results.all),graph.signaling)
#'
#' ## Signaling collection
#'
#' if(interactive()){
#' cytoID.signaling<-sourceCytoscape(results.all,
#'     name.graphs = graph.signaling, collection.name ="SignalingPathway")
#' }
#'
#' ## Other collection
#' if(interactive()){
#' cytoID.other<-sourceCytoscape(results.all,
#'     name.graphs = graph.other, collection.name ="OtherPathway")
#' }
#'
#'
#' @export
sourceCytoscape<-function(sourceObj,
                          name.graphs=names(sourceObj),
                          collection.name="SourceCollection",
                          map.name.variable=NULL,
                          method="bonferroni"){

  ### Gli indici sono calcolati su tutti i grafi
  # aggiungi warnings:
  # per calcolarli su un sottoinsieme, ripetere l'analisi

  ### Verify connection of Cytoscape
  connection<-tryCatch(r2cytoscape::checkCytoscapeVersion(),
                       error=function(c){
                         stop("Cytoscape is not connected: lunch cytoscope before call sourceCytoscape function")
                         return(0)
                       })
  if(length(connection)>=2) message(paste0("Cytoscape ",connection[[2]]," (apiVersion ",connection[[1]],") connected"))

  ### check sourceObj type
  # da fare

  ### Check name.graphs
  ind.g<-which(name.graphs %in% names(sourceObj))
  ind.not<-setdiff(1:length(name.graphs),ind.g)
  if(length(ind.not)>0){
    msg<-paste0("One or more graph names are not in sourceObj (i.e., ",
                paste(name.graphs[ind.not],collapse = ", "),"): check the spelling")
  }
  sO<-sourceObj[name.graphs[ind.g]]

  if(length(sO)==0) stop(msg)
  else if(length(ind.not)>0) warning(msg)

  ### Compute node statistics
  cat(paste0("Compute node statistics (based on ",length(sourceObj)," analyzed graphs)")) # based on...
  info.variable<-infoSource(sourceObj = sourceObj,method= method)$variable

  ind.Inf<-which(info.variable$score==Inf)
  if(length(ind.Inf)>0) info.variable$score[ind.Inf]<-.Machine$double.xmax
  cat("\n\n")

  ### Check map.name.variable (se non riesco a mappare un nome utilizzo di default 'none')
  if(!is.null(map.name.variable)){
    ind.name.var<-match(rownames(info.variable),names(map.name.variable))
    n.na<-sum(is.na(ind.name.var))
    if(n.na>0){
      warn.msg<-paste0("One or more variable names are not mapped in map.name.variable (n=",n.na,"): \'none\' name is used")
      warning(warn.msg)
      mpn<-append(map.name.variable,list("none"))
      names(mpn)[length(mpn)]<-"none"
      ind.name.var[is.na(ind.name.var)]<-length(mpn)
    } else { mpn<-map.name.variable }
    mapped.name<-mpn[ind.name.var]
    mapped.name.vec<-as.vector(sapply(mapped.name,paste,collapse=";"))
  } else { mapped.name.vec<- rownames(info.variable) }
  info.variable<-data.frame("mapped.name"=mapped.name.vec,info.variable,stringsAsFactors = FALSE)


  #### STEP 1: load graphs and data tables in cytoscape

  # network ids
  network.ids<-rep(NA,length(sO))
  network.names<-rep(NA,length(sO))

  for(i in 1:length(sO)){

    # convert graph
    ig<-igraph::igraph.from.graphNEL(sO[[i]]$Graph)

    # edges table
    tab.edges<-igraph::as_data_frame(x = ig,what = "edges")
    colnames(tab.edges)<-c("source","target","weight")

    # nodes table
    tab.nodes<-igraph::as_data_frame(x=ig,what = "vertices")
    source<-rep(0,nrow(tab.nodes))
    source[tab.nodes$name %in% unique(unlist(sO[[i]]$orderingSet))]<-1
    source[tab.nodes$name %in% sO[[i]]$primarySet]<-2

    # info nodes table
    ind.info<-match(tab.nodes$name,rownames(info.variable))
    tab.nodes<-data.frame("id"=tab.nodes$name,
                          source=source,
                          info.variable[ind.info,],
                          #node.size,border.paint,font.size,
                          row.names = NULL,stringsAsFactors = FALSE)


    ## Remove special character
    net.name<-iconv(names(sO)[i], to = "ASCII//TRANSLIT")
    net.name <- gsub(" ","",tools::toTitleCase(net.name))

    ## create network
    cat(paste0(">>> Load: ",net.name," <<<\n"))
    network.ids[i] <- r2cytoscape::createNetwork(edges = tab.edges,network.name = net.name,collection.name = collection.name)
    network.names[i]<-net.name
    ## load nodes table
    r2cytoscape::loadTableData(tab.nodes,data.key.column = "id")
    cat("\n")

  }

  #### STEP 2: import default source set style
  path.style<-system.file("extdata","SourceSetStyle.xml",package = "SourceSet")
  style.name<-r2cytoscape::commandRun(paste0("vizmap load file file=",path.style))
  r2cytoscape::updateStyleMapping(style.name[2], r2cytoscape::mapVisualProperty("node fill color","relevance","c",c(0,max(info.variable$relevance)),c("#F7F7F7","#2C60FA")))



  #### STEP 3: apply SourceSet style to each graphs
  cat("Update style..\n")
  msg<-sapply(network.ids,function(id){
    r2cytoscape::setCurrentNetwork(network = id)
    Sys.sleep(0.1)
    r2cytoscape::applyStyle(style.name[2])
  })

  return(list(
    "collection"=collection.name,
    "network"=network.names,
    "cytoscape.id"=network.ids
  ))

}

#' Visualize in Cytoscape the graphical union induced by the source sets of a collection of graphs
#'
#' The function, thanks to the connection with the Cytoscape software, allows the user to create the graphical union induced by the source sets of a collection of graphs to be visualized in a unique session, while documenting interesting findings.
#' @param sourceObj a \code{SourceSetObj} objects, i.e. the output of the \code{\link[SourceSet]{sourceSet}} function.
#' @param name.graphs the names of the graphs to be visualized. Default value is \code{names(sourceObj)}. NOTE: even if a subset of graphs are selected in \code{name.graphs}, the returned statistics are always calculated on the entire collection in the \code{sourceObj} argument.
#' @param collection.name name of the collection of graphs displayed in Cytoscape.
#' @param network.name name of the resulting union graph.
#' @param map.name.variable a list of customized labels to be associated with the names of the genes. Each list element must contain only one value (i.e. the new label), and the name of each element must be associated with the names of the genes given as input to the \code{\link[SourceSet]{sourceSet}} function (column names of \code{data} input argument). If a label is not mapped, the original name is used.
#' @param method correction method for p-values calculated on graphs. The adjustment methods allowed are: \code{fdr} (default), \code{holm}, \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH}, \code{BY} or \code{none}. For more details refer to \code{\link[stats]{p.adjust}}.
#' @param complete.edges if \code{TRUE}, the graphs selected in \code{name.graphs} are merged and the induced graph of the variables that appear at least in a source set is returned. if \code{FALSE}, the subgraph induced by the variables in the source set of each graph specified in \code{name.graphs} is found, and their union is returned.
#' @param return.unionGraph if \code{TRUE}, the function returns the data frame of the edges of the resulting union graph, together with information about the variables obtained internally through the \code{\link[SourceSet]{infoSource}} function.
#' @details The visual node attributes size and fill color are defined in a dynamic manner through a visual mapping based on the indices provided by the \code{\link[SourceSet]{infoSource}} function (automatically uploaded in the bottom panel - right side).
#'
#' A continous mapper between \code{sub.n.source} attribute and size is applied: higher values are represented with bigger nodes.
#' On the other hand, a color gradient mapper between fill node color and \code{relevance} is adopted: higher values are highlighted with darker blue color.
#'
#' The edges connecting nodes belonging to the graph induced by the source set of each graph are represented by a solid line; while, the edges that connect two variables linked in the union of the graphs, but not within the same source set of a single graph, have dotted lines (supplied only if \code{complete.edges=TRUE}).
#'
#' The default style can be changed manually either within Cytoscape (for further information see \href{http://manual.cytoscape.org/en/stable/Styles.html}{manual}) or within an R package \code{r2cytoscape} through network SUID returned by the \code{sourceCytoscape} function (for further details see \href{https://github.com/cytoscape/r2cytoscape}{manual}).
#'
#' It is also possible to call the sourceCytoscape function multiple times, with all the graphs being visualized in a unique session within a collection specified by collection.name.

#' @note The function use the \code{r2cytoscape} package to connect to Cytoscape from R using CyREST. \code{r2cytoscape} can be downloaded from:
#' \itemize{
#'    \item{Bioconductor: \code{biocLite("r2cytoscape")};}
#'    \item{GitHub: \code{install_github("cytoscape/r2cytoscape")}.}
#' }
#'
#' To enable the display function to work properly, three simple steps are required:
#' \itemize{
#'    \item{ Download \href{https://cytoscape.org/download.php}{Cytoscape} (version 3.3 or later);}
#'    \item{ Complete installation wizard;}
#'    \item{ Launch Cytoscape (before calling the functions).}
#' }
#' @seealso \code{\link[SourceSet]{sourceSet}}, \code{\link[SourceSet]{sourceCytoscape}}, \code{r2cytoscape}
#' @examples
#' ## Load the SourceSetObj obtained from the source set analysis of ALL dataset
#'
#' # see vignette for more details
#' print(load(file=system.file("extdata","ALLsourceresult.RData",package = "SourceSet")))
#' class(results.all)
#'
#' ## NB: Remember to launch cytoscape before running the following commands
#' # Create two collections of pathways to visualize the results
#' graph.signaling<-names(results.all)[grep("signaling",names(results.all))]
#' graph.other<-setdiff(names(results.all),graph.signaling)
#'
#' ## Signaling collection
#'
#' if(interactive()){
#' cytoID.signaling.union<-sourceUnionCytoscape(results.all,
#'       name.graphs =graph.signaling ,collection.name ="SignalingPathway",
#'       network.name ="SignalingUnion")
#' }
#' ## Other collection
#'
#' if(interactive()){
#' cytoID.other.union<-sourceUnionCytoscape(results.all ,
#'       name.graphs =graph.other,collection.name ="OtherPathway" ,
#'       network.name ="OtherUnion")
#' }
#' @return The function returns an interactive session in Cytoscape.
#' @export
sourceUnionCytoscape<-function(sourceObj,name.graphs=names(sourceObj),
                                  collection.name="SourceSetUnion",
                                  network.name="UnionSourceSetsGraph",
                                  map.name.variable=NULL,
                                  method="bonferroni",
                                  complete.edges=TRUE,
                                  return.unionGraph=FALSE){

  ### Verify connection of Cytoscape
  connection<-tryCatch(r2cytoscape::checkCytoscapeVersion(),
                       error=function(c){
                         stop("Cytoscape is not connected: lunch cytoscope before call sourceCytoscape function")
                         return(0)
                       })
  if(length(connection)>=2) message(paste0("Cytoscape ",connection[[2]]," (apiVersion ",connection[[1]],") connected"))


  ### Check name.graphs
  ind.g<-which(name.graphs %in% names(sourceObj))
  ind.not<-setdiff(1:length(name.graphs),ind.g)
  if(length(ind.not)>0){
    msg<-paste0("One or more graph names are not in sourceObj (i.e., ",
                paste(name.graphs[ind.not],collapse = ", "),"): check the spelling")
  }
  sO<-sourceObj[name.graphs[ind.g]]

  if(length(sO)==0) stop(msg)
  else if(length(ind.not)>0) warning(msg)


  ### Check if exist at least one non empty source set
  ind.empty<- which(sapply(sO,function(x) length(x$primarySet))==0)
  if(length(ind.empty)>0) sO<-sO[-ind.empty]
  if(length(sO)==0) stop("The source set is empty for all the analyzed graphs")

  sub.union.source<-table(unlist(lapply(sO,function(x) x$primarySet)))
  union.source<-names(sub.union.source)
  n.graph<-length(sO)

  ### Compute node statistics
  cat(paste0("Compute node statistics (based on ",length(sourceObj)," analyzed graphs)")) # based on...
  info.variable<-infoSource(sourceObj = sourceObj,method=method)$variable
  ind.Inf<-which(info.variable$score==Inf)
  if(length(ind.Inf)>0) info.variable$score[ind.Inf]<-.Machine$double.xmax
  cat("\n\n")


  ### Check map.name.variable (se non riesco a mappare un nome utilizzo di default 'none')
  if(!is.null(map.name.variable)){
    ind.name.var<-match(rownames(info.variable),names(map.name.variable))
    n.na<-sum(is.na(ind.name.var))
    if(n.na>0){
      warn.msg<-paste0("One or more variable names are not mapped in map.name.variable (n=",n.na,"): \'none\' name is used")
      warning(warn.msg)
      mpn<-append(map.name.variable,list("none"))
      names(mpn)[length(mpn)]<-"none"
      ind.name.var[is.na(ind.name.var)]<-length(mpn)
    } else { mpn<-map.name.variable }
    mapped.name<-mpn[ind.name.var]
    mapped.name.vec<-as.vector(sapply(mapped.name,paste,collapse=";"))
  } else { mapped.name.vec<- rownames(info.variable) }
  info.variable<-data.frame("mapped.name"=mapped.name.vec,info.variable,stringsAsFactors = FALSE)


  ## check union param
  #if( !(union %in% c("partial","complete")) ) stop("union parameter must be \'partial\' or \'complete\'")

  M_source<-M_ext<-matrix(0,nrow = length(union.source), ncol = length(union.source),
                          dimnames = list(union.source,union.source))


  for(i in 1:n.graph){
    mi<-methods::as(sO[[i]]$Graph,"matrix")

    ### edges between nodes that belong to the source set of:
    ## the graph i
    mi_source<-mi[sO[[i]]$primarySet,sO[[i]]$primarySet]

    ## other graph
    ind.mi<-which(rownames(mi) %in% union.source)
    mi_ext<-mi[ind.mi,ind.mi]


    M_source[rownames(mi_source),colnames(mi_source)]<-M_source[rownames(mi_source),colnames(mi_source)]+mi_source
    M_ext[rownames(mi_ext),colnames(mi_ext)]<-M_ext[rownames(mi_ext),colnames(mi_ext)]+mi_ext
  }

  ### Edge tables
  # df.source: 1. induced graph by source set 2.union of source sets
  # df.ext: 1. union of the graphs 2. induced graph by source sets
  ii<-which(upper.tri(M_ext),arr.ind = TRUE)
  df.source<-data.frame("source"=union.source[ii[,1]],"target"=union.source[ii[,2]],"weight"=M_source[ii],stringsAsFactors = FALSE)
  df.ext<-data.frame("source"=union.source[ii[,1]],"target"=union.source[ii[,2]],"weight"=M_ext[ii],stringsAsFactors = FALSE)

  ## complete=TRUE: plot both type of edges; complete=FALSE: only proper edges
  if(complete.edges) ind.no.edge<-which(df.ext$weight==0)
  else ind.no.edge<-which(df.source$weight==0)

  if(length(ind.no.edge)>0){
    df.source<-df.source[-ind.no.edge,]
    df.ext<-df.ext[-ind.no.edge,]
  }


  ## type edge: solid (1.proper source set edge); dashed (0.nodes in union source set connected in at least one graph)
  edge.type<-as.numeric(df.source$weight>0)

  if(complete.edges) tab.edges<-data.frame(df.ext,"type"=edge.type,stringsAsFactors = FALSE)
  else tab.edges<-data.frame(df.source,"type"=edge.type,stringsAsFactors = FALSE)


  # info nodes table
  ind.info<-match(union.source,rownames(info.variable))
  tab.nodes<-data.frame("id"=union.source,
                        sub.n.source=as.numeric(sub.union.source),
                        info.variable[ind.info,],
                        #node.size,border.paint,font.size,
                        row.names = NULL,stringsAsFactors = FALSE)

  id.nodes<-tab.nodes[,"id",drop=FALSE]


  ## create network
  cat(paste0(">>> Load: ",network.name," <<<\n"))
  id <- r2cytoscape::createNetwork(nodes = id.nodes,edges = tab.edges,network.name = network.name,collection.name = collection.name)
  ## load nodes table
  r2cytoscape::loadTableData(tab.nodes,data.key.column = "id")
  cat("\n")


  #### STEP 2: import default source set style
  path.style<-system.file("extdata","SourceSetUnionStyle.xml",package = "SourceSet")
  style.name<-r2cytoscape::commandRun(paste0("vizmap load file file=",path.style))

  r2cytoscape::updateStyleMapping(style.name[2], r2cytoscape::mapVisualProperty("node fill color","relevance","c",c(0,max(info.variable$relevance)),c("#F7F7F7","#2C60FA")))
  r2cytoscape::updateStyleMapping(style.name[2], r2cytoscape::mapVisualProperty("node size","sub.n.source","c",c(0,max(tab.nodes$sub.n.source)),c(30,100)))
  r2cytoscape::updateStyleMapping(style.name[2], r2cytoscape::mapVisualProperty("node label font size","sub.n.source","c",c(0,max(tab.nodes$sub.n.source)),c(15,30)))

  r2cytoscape::setCurrentNetwork(network = id)
  r2cytoscape::applyStyle(style.name[2])


  res<-list(
    "collection"=collection.name,
    "network"=network.name,
    "cytoscape.id"=id
  )


  if(return.unionGraph) {
    row.names(tab.edges)<-NULL
    unionGraph<-list(list("edges"=tab.edges,"nodes"=tab.nodes))
    names(unionGraph)<-"unionGraph"
    res<-append(res,unionGraph)
  }

  return(res)
}






