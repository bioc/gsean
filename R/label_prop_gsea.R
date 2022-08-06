# label propagation for GSEA
label_prop_gsea <- function(geneset, x, adjacency, threshold = 0.99, nperm = 1000,
                            minSize = 1, maxSize = Inf, gseaParam = 1, nproc = 0,
                            BPPARAM = NULL, ...)
{
  Rownames <- rownames(adjacency)
  Colnames <- colnames(adjacency)
  if(!is.null(Rownames) & !is.null(Colnames))
    if(!identical(Rownames, Colnames))
      stop("the row or column names of the adjacency matrix are not identical")
  if(is.null(Rownames) & is.null(Colnames))
    stop("the row or column names of the adjacency matrix are required")
  
  Names <- Rownames
  if(is.null(Names)) Names <- Colnames
  x <- na.omit(match(x, Names))
  if(length(x) == 0)
    stop("No genes in both the list and the adjacency matrix")
  
  adjacency[adjacency >= threshold] = 1
  adjacency[adjacency < threshold] = 0
  lp <- label.prop(adjacency, ind.positives = x, ...)
  scaled.scores <- as.vector(scale(lp$p))
  names(scaled.scores) <- Names
  statistic <- scaled.scores
  
  result.GSEA <- fgsea(geneset, statistic, nperm = nperm, minSize = minSize,
                       maxSize = maxSize, nproc = nproc, gseaParam = gseaParam,
                       BPPARAM = BPPARAM)
  result.GSEA
}


# from the package RANKS

#' Function that implements the Label propagation algorithm of Zhu and Ghahramani
#' @param W : adjacency matrix of the graph
#' @param ind.positives : indices of the "core" positive examples of the graph.
#'                They represent to the indices of W corresponding to the positive examples
#' @param tmax : maximum number of iterations (def: 1000)
#' @param eps : maximum allowed difference between the computed probabilities at the steady state (def. 1e-5)
#' @param norm : if TRUE (def) the adjacency matrix W of the graph is normalized to M = D^-1 * W, otherwise
#'        it is assumed that the matrix W is just normalized
#' @return a list with three elements:
#' - p : the probability at the steady state
#' - ind.positives: indices of the "core" positive examples of the graph (it is equal to the same
#'                  input parameter)
#' - n.iter : number of performed iterations
label.prop <- function(W, ind.positives, tmax=1000, eps=1e-5, norm=TRUE) { 
  if (norm) 
    M <-Prob.norm(W) # M = D^-1 * W
  else
    M <- W;
  n <- nrow(M);
  p <- numeric(n);
  names(p)  <- rownames(W);
  n.positives <- length(ind.positives);
  if (n.positives == 0)
    stop("label.prop: number of core positives is equal to 0!");
  p[ind.positives] <- 1;
  
  M <- t(M);
  for (t in 1:tmax) {
    pold <- p;
    p <- M %*% pold;
    p[ind.positives] <- 1;
    if (norm1(p-pold) < eps) break();  
  }  
  return(list(p=p, ind.positives=ind.positives, n.iter=t));  
}


# from the package NetPreProc

########################################################
setGeneric("Prob.norm", 
           function(W) standardGeneric("Prob.norm"));

# Probabilistic normalization of a graph.
# Method to normalize weights of square symmetric adjacency matrices
# A network matrix is normalized by dividing each entry W_ij by the the sum of elements of row i 
# In other words if D is a diagonal matrix such that D_ii = \sum_j W_ij then the normalize matrix is: W_norm = D^(-1) * W 
# Input:
# W : symmetric matrix 
# Output:
# a normalized  matrix
# N.B.: La matrice risultante non e' simmetrica!
setMethod("Prob.norm", signature(W="matrix"),
          function(W) {
            n <- nrow(W);
            if (n != ncol(W))
              stop("first arg must be a square matrix");
            names.var <- rownames(W);
            diag.D <- apply(W,1,sum);
            diag.D[diag.D==0] <- Inf;
            inv.diag.D <- 1/diag.D;
            W <- .C("norm_2", as.double(W), as.double(inv.diag.D), as.integer(n), PACKAGE="NetPreProc")[[1]];
            W <- matrix(W, nrow=n);
            rownames(W) <- colnames(W) <- names.var;
            return(W);
          })

# Probabilistic normalization of a graph.
# Method to normalize weights of square symmetric adjacency matrices
# A network matrix is normalized by dividing each entry W_ij by the the sum of elements of row i 
# In other words if D is a diagonal matrix such that D_ii = \sum_j W_ij then the normalize matrix is: W_norm = D^(-1) * W 
# Input:
# W : an object of the virtual class graph (hence including objects of class graphAM  and graphNEL from the package graph).
# Output:
# a normalized  matrix
# N.B.: La matrice risultante non e' simmetrica!
setMethod("Prob.norm", signature(W="graph"),
          function(W) {
            W <- as(W, "matrix");
            return(Prob.norm(W));
          })
