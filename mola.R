
# set seed for reproducibility
set.seed(123)

# set up packages
if(!require("reticulate")) install.packages("reticulate")
if(!require("dplyr")) install.packages("dplyr")
if(!require("vegan")) install.packages("vegan")
if(!require("nlme")) install.packages("nlme")
if(!require("survival")) install.packages("survival")
library(reticulate)
library(dplyr)
library(vegan)
library(nlme)
library(survival)

# create or load environment
envname <- "matilda-env"

if (!virtualenv_exists(envname)) {
  virtualenv_create(envname = envname, python = "3.13")
  use_virtualenv(envname, required = TRUE)
  
  # install packages
  py_install(
    packages = "../matilda", # do ../matilda-main if downloaded zip file
    pip = TRUE,
    pip_options = c("--upgrade", "--no-user")
  )
  py_install(
    packages = "pandas",
    pip = TRUE,
    pip_options = c("--upgrade", "--no-user")
  )
  
}else
{
  # just activate environment
  use_virtualenv(envname, required = TRUE)
  # py_config() # to verify 
}

# import python modules
matilda <- import("matilda")
harmonic <- import("matilda.harmonic", convert = FALSE)

py_run_string("
from collections.abc import Iterable
import numpy as np

import numpy as np
import pandas as _pd
import math as _math
import re as _re

from numpy import ones, ndarray, percentile, array, float64, log, exp, where, concatenate, reshape, delete, inf, intersect1d, apply_along_axis, append, empty
from numpy.testing import assert_almost_equal
from itertools import repeat, combinations
from random import sample, choices
from inspect import getfullargspec
from scipy.special import digamma
import warnings
from scipy.spatial.distance import cdist

# functions to check persistence diagrams and convert diagrams between formats
def check_diagram(D):
    '''Checks for persistence diagrams.
    
    Internal method to verify that birth values are non-negative and less than death values.

    Parameters
    ----------
    `D` : numpy.ndarray
        The input diagram to be checked.
    
    Returns
    -------
    None
    '''
    # make sure that birth values are <= death values
    # if D has one entry
    if len(D.shape) != 1:
        pers = D[:,1] - D[:,0]
        birth = D[:,0]
    else:
        pers = D[1] - D[0]
        birth = D[0]
    if where(pers < 0)[0].shape[0] > 0:
        raise Exception('Birth values have to be <= death values.')
    # make sure that birth values are >= 0
    if where(birth < 0)[0].shape[0] > 0:
        raise Exception('Birth values have to be non-negative.')

def preprocess_diagram(D, inf_replace_val = None, ret = False):
    '''Verify the format of a persistence diagram and convert to a standard format.
    
    This function can verify a persistence diagram from the ripser, gph, flagser, gudhi or cechmate packages 
    and convert any such diagram into a list of numpy arrays if desired (largely an internal functionality but
    can be used in a standalone fashion).
    
    Parameters
    ----------
    `D` : any
        The persistence diagram to be verified. An exception will be raised if `D` is not a persistence
        diagram computed from one of the aforementioned packages.
    `inf_replace_val` : float or int, default `None`
        The value with which `inf` values should be replaced, if desired.
    `ret` : bool, default `False`
        Whether or not to return a processed diagram.
    
    Returns
    -------
    None or list of numpy.ndarray
        If `ret` is `True` and the diagram is verified then a list is returned. The i-th element of 
        the returned list is the array of i-dimensional topological features in the diagram.
    '''
    # error check inf_replace_val
    def check_val(D, inf_replace_val):
        if inf_replace_val != None:
            if not isinstance(inf_replace_val, type(2.0)) and not isinstance(inf_replace_val, type(2)):
                raise Exception('inf_replace_val must be a number.')
            if inf_replace_val <= 0:
                raise Exception('inf_replace_val must be positive.')
            max_death = max([d[d[:,1] != float('inf')].max() if len(d[d[:,1] != float('inf')]) > 0 else 0 for d in D])
            if inf_replace_val < max_death:
                raise Exception('inf_replace_val should be at least as large as any death value in the diagram.')
            # replace vals
            D[0][D[0][:,1] == float('inf'),1] = inf_replace_val
            return D
        else:
            return D
    error_message = 'Diagrams must be computed from either the ripser, gph, flagser, gudhi or cechmat libraries.'
    # first check if the diagram is from ripser, gph or flagser
    if isinstance(D, dict):
        if not 'dgms' in D.keys():
            raise Exception(error_message)
        if not isinstance(D['dgms'],list):
            raise Exception(error_message)
        if set([type(x) for x in D['dgms']]) != set([type(array([0,1]))]):
            raise Exception(error_message)
        if set([len(x.shape) for x in D['dgms']]) != set([2]):
            raise Exception(error_message)
        if set([x.shape[1] for x in D['dgms']]) != set([2]):
            raise Exception(error_message)
        # convert to list of numpy arrays
        D = D['dgms']
        # perform final numeric diagram checks
        lst = [check_diagram(d) for d in D]
        if ret == True:
            return check_val(D, inf_replace_val)
    # now check if diagram is from cechmate or gudhi
    elif isinstance(D, list):
        if set([type(x) for x in D]) == set([type((1,2))]):
            if set([len(x) for x in D]) != set([2]):
                raise Exception(error_message)
            if set([type(x[1]) for x in D]) != set([type((1,2))]):
                raise Exception(error_message)
            if set([len(x[1]) for x in D]) != set([2]):
                raise Exception(error_message)
            dims = [x[0] for x in D]
            if any(e < 0 for e in dims):
                raise Exception(error_message)
            if set([type(d) for d in dims]) != set([type(1)]):
                raise Exception(error_message)
            # convert to list of numpy arrays
            max_dim = max(dims)
            res = [array([0,0]).reshape((1,2)) for i in range(max_dim + 1)]
            for pt in D:
                res[pt[0]] = append(res[pt[0]], array([pt[1][0], pt[1][1]]).reshape((1,2)),axis = 0)
            res = [r[range(1,len(r)),:] for r in res]
            lst = [check_diagram(d) for d in res]
            if ret == True:
                return check_val(res, inf_replace_val)
        elif set([type(x) for x in D]) == set([type(array([0,1]))]):
            if set([len(x.shape) for x in D]) != set([2]):
                raise Exception(error_message)
            if set([x.shape[1] for x in D]) != set([2]):
                raise Exception(error_message)
            # no need to convert as this is the base format
            # final numeric diagram checks
            lst = [check_diagram(d) for d in D]
            if ret == True:
                return check_val(D, inf_replace_val)
        else:
            raise Exception(error_message)
    else:
        raise Exception(error_message)
    
def fdr(pvals, alpha: float):
    pvals = np.asarray(pvals)
    m = pvals.size
    if m == 0:
        return np.zeros(0, dtype=bool)
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    thresholds = alpha * (np.arange(1, m + 1) / m)
    is_sig_sorted = sorted_p <= thresholds
    if not is_sig_sorted.any():
        return 0.0
    k = np.max(np.where(is_sig_sorted))
    cutoff = sorted_p[k]
    return cutoff

class universal_null:
    def __init__(self, dims:list = [1], distance_mat:bool = False, alpha:float = 0.05, infinite_cycle_inference:bool = False):

        if not isinstance(dims, type([0,1])):
            raise Exception('dims must be a list.')
        if set([type(d) for d in dims]) != set([type(1)]):
            raise Exception('Each dimension in dims must be an integer.')
        if min(dims) < 1:
            raise Exception('Each dimension in dims must be at least 1.')
        self.dims = dims

        if isinstance(distance_mat, type(True)) == False:
            raise Exception('distance_mat must be True or False.')
        self.distance_mat = distance_mat

        if not isinstance(alpha,type(0.05)):
            raise Exception('alpha must be a float.')
        if alpha <= 0 or alpha > 1:
            raise Exception('alpha must be between 0 and 1.')
        self.alpha = alpha

        if isinstance(infinite_cycle_inference, type(True)) == False:
            raise Exception('infinite_cycle_inference must be True or False.')
        self.infinite_cycle_inference = infinite_cycle_inference
    def __str__(self):
        '''Describe a universal null procedure based on the dimensions being analyzed,
        whether or not the input will be a distance matrix, the Type 1 error rate (alpha) and
        whether or not infinite cycle inference will be carried out.'''
        dms = ''
        if not self.distance_mat:
            dms = 'non-'
        if len(self.dims) == 1:
            dim_str = 'dimension ' + str(self.dims[0])
        else:
            dim_str = ', '.join([str(d) for d in self.dims[::-1]]) + ' and ' + str(self.dims[-1])
        if self.distance_mat:
            distmat_str = 'distance-matrix'
        else:
            distmat_str = 'point-cloud'
        if self.infinite_cycle_inference:
            ici_str = ''
        else:
            ici_str = 'no '
        s = f'Universal null procedure for {dim_str}, {distmat_str} input, a Type 1 error rate of {str(self.alpha)} and {ici_str}infinite cycle inference.'
        return s
    def compute(self, homology_computer, thresh):
        # create dim 1 persistence diagram
        bars_dim = homology_computer.bars.get(1, {})
        diagram = [concatenate([reshape(array(x), (1, 2)) for x in bars_dim.values()], axis = 0)]
        keys = list(bars_dim.keys())
        # error check the persistence diagram
        try:
            diagram = preprocess_diagram(D = diagram, ret = True)
        except Exception as ex:
            raise Exception('The output of diagam_fun(X, thresh) was not in the correct format for a persistence diagram.')
        # set up return dict
        ret = {}
        # check if there is any work to be done
        if self.alpha == 1:
            return {'subsetted_keys':keys, 'p_values':[1 for _ in range(len(keys))], 'thresh':1.0, 'diagram':diagram, 'keys':keys}
        if len(diagram) > 0:
            # subset for diagrams above dimension 0
            diag_highdim = diagram # just in our H1 case
            # replace inf values with thresh
            diag_highdim = [where(d == inf,concatenate([reshape(d[:,0], (len(d), 1)), reshape(array([thresh for _ in range(len(d))]), (len(d), 1))], axis = 1),d) for d in diag_highdim]
            # compute test statistics and p-values
            A = 1 # for VR persistent homology
            lambd = -1*digamma(1)
            pi = [diag_sub[:,1]/diag_sub[:,0] for diag_sub in diag_highdim] # ratio of death to birth
            loglog_pi = [log(log(pi_sub)) for pi_sub in pi]
            Lbar = [loglog_pi_sub.mean() for loglog_pi_sub in loglog_pi]
            B = -1*lambd - A*Lbar
            test_statistics = [A*loglog_pi[i] + B[i] for i in range(len(loglog_pi))]
            p_values = [exp(-1*exp(test_statistics_sub)) for test_statistics_sub in test_statistics]
            # determine FDR thresholds in each dimension
            alpha_thresh = fdr(p_values[0],self.alpha) if len(p_values[0]) > 0 else Inf
            # subset the keys and return
            subsetted_keys = [keys[i] for i in range(len(keys)) if p_values[0][i] <= self.alpha]
            return {'subsetted_keys':subsetted_keys, 'p_values':p_values[0], 'thresh':alpha_thresh, 'diagram':diag_highdim[0], 'keys':keys}

class FilteredSimplicialComplex(object):

    def __init__(
        self, dimension=0, simplices=None, simplices_indices=None, appears_at=None
    ):
        self.dimension = dimension
        if simplices_indices is None:
            self.simplices_indices = []
        else:
            self.simplices_indices = simplices_indices
        if appears_at is None:
            self.appears_at = []
        else:
            self.appears_at = appears_at
        if simplices is None:
            self.simplices = []
        else:
            self.simplices = simplices

    def get_boundary_matrix_at_filtration(
        self,
        boundary_dict,
        value=None,
        simplex_id=None,
        dim=1,
        mode='economic',
        cofaces_to_skip=[],
        verbose=False,
    ):
        if value:
            indices_list = self.subcomplex_at_filtration(value)
        elif simplex_id:
            indices_list = self.subcomplex_at_index(simplex_id)
        else:
            raise ValueError('specify either value or simplex_id')
            
        # new code
        valid_indices = [
            i for i in indices_list
            if isinstance(self.simplices[self.simplices_indices[i]], Iterable)
            and not isinstance(self.simplices[self.simplices_indices[i]], (int, np.integer))
        ]

        if mode == 'economic':
            faces_idx = np.array(
                [
                    i
                    for i in valid_indices
                    if len(self.simplices[self.simplices_indices[i]]) == dim
                ]
            )
        elif mode == 'full':
            # returns all (d-1)-dimensional simplices ids
            faces_idx = np.array(
                [
                    i
                    for i in valid_indices
                    if len(self.simplices[self.simplices_indices[i]]) == dim
                ]
            )
        else:
            raise ValueError('unknown mode')

        cofaces_idx = np.array(
            [
                i
                for i in valid_indices
                if (len(self.simplices[self.simplices_indices[i]]) == dim + 1)
                and (i not in cofaces_to_skip)
                and (i in boundary_dict)
            ]
        )

        nrows = len(faces_idx)
        ncols = len(cofaces_idx)
        if verbose:
            print('creating {}x{} boundary matrix'.format(nrows, ncols))

        boundary_matrix = np.zeros((nrows, ncols), dtype=int)

        for j, cid in enumerate(cofaces_idx):
            for fid in boundary_dict[cid]:
                i = np.where(faces_idx == fid)[0]
                boundary_matrix[i, j] = boundary_dict[cid][fid]

        return boundary_matrix, faces_idx, cofaces_idx

    def subcomplex_at_filtration(self, t):
        indices_list = []

        for i, id in enumerate(self.simplices_indices):
            if self.appears_at[id] <= t:
                indices_list.append(i)

        return indices_list

    def subcomplex_at_index(self, id):
        import warnings

        if id < 0:
            warnings.warn('{} is less than 0 - returning all simplices')
            return [i for i in self.simplices_indices]
        else:
            return [i for i in self.simplices_indices if i <= id]
            
def _vertex_weights_from_edge_chain(K, edge_chain, n_vertices):
    w = np.zeros(int(n_vertices), dtype=float)
    items = edge_chain.items() if hasattr(edge_chain, 'items') else edge_chain
    max_coef = 0
    for k in edge_chain.keys():
        a = abs(float(edge_chain[k]))
        i = K.simplices[int(k)][0]
        j = K.simplices[int(k)][1]
        if a > max_coef:
            max_coef = a
        w[i] = w[i] + a
        w[j] = w[j] + a

    # normalize (optional; mirrors common notebook practice)
    if max_coef > 0:
        w = w / max_coef
    return w
    
def _get_edge_weights(K, edge_chain):
    ret = np.empty((0, 3))
    for k in edge_chain.keys():
        n = abs(float(edge_chain[k]))
        arr = np.reshape(np.array([K.simplices[int(k)][0], K.simplices[int(k)][1], n]), (1, 3))
        ret = np.concatenate([ret, arr], axis = 0)
    return ret # unnormalized, but edge weights are maximum value 1

              ")

# correlation distance matrix with missing values
get_dist_mat <- function(df){
  
  D <- 1 - cor(t(as.matrix(df)))
  return(D)
  
}

# main wrapper function
run_harmonic_PH <- function(df, metr = 'correlation', alpha=0.2, log_str, log_filepath, top = NULL){
  tryCatch(expr = {
    
    # create return list
    ret_list <- list()
    
    # compute distance matrix
    if(metr == "correlation")
    {
      D <- get_dist_mat(df)
    }else
    {
      D <- as.matrix(dist(df))
    }
    
    log_str <- paste0(log_str, 'Computed distance matrix at ', Sys.time(), ': dimension (', dim(D)[[1]],',',dim(D)[[2]],').\n')
    
    enc_rad <- min(apply(D, MARGIN = 1L, FUN = max))
    
    log_str <- paste0(log_str, 'Computed enclosing radius at ', Sys.time(), '. Value: ', enc_rad,'\n')
    
    ret_list$enc_rad <- enc_rad
    
    # compute homology
    max_dim <- 2L
    K <- matilda$FilteredSimplicialComplex()
    K$construct_vietoris_from_metric(D, max_dim, enc_rad)
    hom <- matilda$PersistentHomologyComputer()
    hom$compute_persistent_homology(K, with_representatives=T, modulus=0L)
    log_str <- paste0(log_str, 'Computed PH at ', Sys.time(), '. Number of H1 bars: ', length(hom$bars$`1`),'\n')
    
    if(length(hom$bars$`1`) == 0)
    {
      ret_list$diagram <- data.frame(birth = numeric(), death = numeric(), p_values = numeric(), keys = numeric())
      ret_list$point_weights_mat <- matrix(data = 0, ncol = 0, nrow = 0)
      ret_list$edge_weights_list <- list()
      ret_list$alpha_thresh = 0.0
      ret_list$p_values = list()
      ret_list$log_str <- paste0(log_str, 'Finished calculations at ', Sys.time())
      return(ret_list)
    }
    
    # filter for FDR
    univ_null <- py$universal_null(distance_mat = T, alpha = alpha)
    univ_null_results_list <- univ_null$compute(hom, enc_rad)
    
    log_str <- paste0(log_str, 'Computed universal null at ', Sys.time(), '. Number of significant features: ', length(univ_null_results_list$subsetted_keys), ' , number of features: ', length(hom$bars$`1`),'\n')
    
    ret_list$diagram <- univ_null_results_list$diagram
    ret_list$keys <- univ_null_results_list$keys
    ret_list$alpha_thresh <- univ_null_results_list$thresh
    ret_list$p_values <- univ_null_results_list$p_values
    
    if(length(univ_null_results_list$subsetted_keys) == 0)
    {
      ret_list$point_weights_mat <- matrix(data = 0,nrow = 0,ncol = 0)
      ret_list$edge_weights_list <- list()
      ret_list$log_str <- paste0(log_str, 'Finished calculations at ', Sys.time())
      return(ret_list)
    }
    
    # copy the contents of K into a new fsc
    K2 <- py$FilteredSimplicialComplex(
      simplices = r_to_py(K$simplices),
      simplices_indices = r_to_py(K$simplices_indices),
      appears_at = r_to_py(K$appears_at)
    )
    
    # Compute harmonic representatives for the selected H1 bars
    if(length(univ_null_results_list$subsetted_keys) == 1)
    {
      univ_null_results_list$subsetted_keys <- list(univ_null_results_list$subsetted_keys)
    }
    
    # subset for top k if desired
    if(!is.null(top))
    {
      diagram <- univ_null_results_list$diagram
      pers <- diagram[, 2L] - diagram[, 1L]
      ordered_pers <- rev(pers[order(pers)])
      thresh <- (ordered_pers[[top]] + ordered_pers[[top + 1]])/2.0
      univ_null_results_list$subsetted_keys <- univ_null_results_list$keys[which(pers > thresh)]
      if(length(which(pers > thresh)) == 1)
      {
        univ_null_results_list$subsetted_keys <- list(univ_null_results_list$subsetted_keys)
      }
      print(univ_null_results_list$subsetted_keys)
    }
    
    log_str <- paste0(log_str, 'Starting harmonic PH at ', Sys.time(), '\n')
    
    hc <- harmonic$HarmonicRepresentativesComputer(K2, hom)
    hc$compute_harmonic_cycles(
      dim = 1L,
      selected_cycles = univ_null_results_list$subsetted_keys,
      precompute = TRUE,
      verbose = 2L,
    )
    
    log_str <- paste0(log_str, 'Computed harmonic PH at ', Sys.time(), '\n')
    
    # get edge and point weight matrices
    n <- nrow(df)
    edge_weights_list <- sapply(
      py_to_r(hc$harmonic_cycles[[1L]]),
      function(edge_chain) py$`_get_edge_weights`(K2, edge_chain)
    )
    point_weights_mat <- sapply(
      py_to_r(hc$harmonic_cycles[[1L]]),
      function(edge_chain) py$`_vertex_weights_from_edge_chain`(K2, edge_chain, n)
    )
    
    # filter for only the desired cycles
    inds <- c()
    for(i in 1:ncol(point_weights_mat)){
      if(as.numeric(colnames(point_weights_mat)[[i]]) %in% univ_null_results_list$subsetted_keys){
        inds <- c(inds, i)
      }
    }
    point_weights_mat <- point_weights_mat[, inds]
    
    if(length(inds) == 1)
    {
      # reformat matrix to have one column of that name
      point_weights_mat <- matrix(data = point_weights_mat, ncol = 1L)
      colnames(point_weights_mat)[[1]] <- as.character(univ_null_results_list$subsetted_keys[[1]])
    }
    
    # also filter edge_weights_list
    inds <- c()
    for(i in 1:length(edge_weights_list)){
      if(length(names(edge_weights_list)) > 0 & length(univ_null_results_list$subsetted_keys) > 0)
      {
        if(as.numeric(names(edge_weights_list)[[i]]) %in% univ_null_results_list$subsetted_keys){
          inds <- c(inds, i)
        }
      }
    }
    if(length(inds) > 0)
    {
      edge_weights_list <- edge_weights_list[inds]
    }else
    {
      edge_weights_list <- list()
    }
    
    ret_list$point_weights_mat <- point_weights_mat
    ret_list$edge_weights_list <- edge_weights_list
    
    ret_list$log_str <- log_str
    
    return(ret_list)
    
  }, error = function(e){
    log_str <- paste0(log_str, 'Error: ',e)
    cat(log_str, file = log_filepath, append = TRUE)
  })
  
}

# run mola on a gene-omic dataframe "df"
run_mola <- function(df, results_dir = '..', alpha = 1.0, verbose = T){
  
  if(verbose)
  {
    print(paste0('Starting MOLA analysis at ', Sys.time(), ' with alpha = ', alpha, '.\n'))
  }
  
  # write to log
  log_str <- paste0('Starting MOLA analysis at ', Sys.time(), ' with alpha = ', alpha, '.\n')
  log_filepath <- paste0(results_dir,'/log.txt')
  cat(log_str, file = log_filepath, append = TRUE)
  
  # run wrapper function
  res <- run_harmonic_PH(df, 'correlation', alpha, log_str, log_filepath)
  
  # check and save features
  point_weights_mat <- res$point_weights_mat
  if(!is.null(point_weights_mat))
  {
    if(ncol(point_weights_mat) > 0)
    {
      save(point_weights_mat, file = paste0(results_dir,'/gene_weights.RData'))
    }
  }else
  {
    if(verbose)
    {
      print("No significant features found.") 
    }
  }
  
  edge_weights_list <- res$edge_weights_list
  if(!is.null(edge_weights_list))
  {
    if(length(edge_weights_list) > 0)
    {
      save(edge_weights_list, file = paste0(results_dir,'/gene_interaction_weights.RData'))
    }
  }
  
  diagram <- data.frame(res$diagram)
  if(ncol(diagram) == 2)
  {
    colnames(diagram) <- c('birth', 'death')
    diagram$p_value <- res$p_values
    diagram$key <- res$keys
    write.csv(diagram, paste0(results_dir,'/persistence_diagram.csv'))
    if(verbose)
    {
      print(paste0('Minimum feature p-value: ', min(res$p_values)))
    }
  }
  
  enc_rad <- res$enc_rad
  alpha_thresh <- res$alpha_thresh
  params_list <- list('enc_rad' = enc_rad, 'alpha_thresh' = alpha_thresh)
  save(params_list, file = paste0(results_dir,'/params.RData'))
  
  cat(res$log_str, file = log_filepath, append = TRUE)
  
}

# filter a (joint sample) diagram for significant loops based on the
# universal null procedure (with a fixed threshold)
filter_loops <- function(diagram, alpha = 0.05){
  
  # define constants
  A <- 1
  lambd <- digamma(1)
  
  # run procedure
  pi <- diagram$death/diagram$birth
  loglog_pi <- log(log(pi))
  Lbar <- mean(loglog_pi, na.rm = T)
  B <- lambd - A*Lbar # negative value
  test_statistics <- A*loglog_pi + B
  pvals <- exp(-1*exp(test_statistics))
  
  # filter
  return(diagram[which(pvals < alpha), ])
  
}

# run a covariate shift analysis on mola results
# "weights" is a matrix with one row per loop, columns being the genes and entries being the weights (from gene_weights.RData files output from run_mola)
# "pheno" is a dataframe with one row per loop containing the sample-level statistics
# "genes" is the vector of shared genes, corresponding to the rows of weights
# "results_dir" is the path of where to save the result file, tvals.RData
covariate_shift_analysis <- function(weights, pheno, genes, results_dir){
  
  # run modelling in a loop
  get_gene_t_for_term <- function(models, term, gene_ids = NULL) {
    tvals <- vapply(models, function(z) {
      ind <- which(names(z) == term)
      return(z[[ind]])
    }, numeric(1))
    if (!is.null(gene_ids)) names(tvals) <- gene_ids
    tvals <- tvals[is.finite(tvals)]
    sort(tvals, decreasing = TRUE)
  }
  
  # fit models predicting the harmonic weights for all genes based on sample-level statistics
  # MODIFY MODELLING FOR YOUR OWN APPLICATION
  dat <- pheno
  models <- lapply(X = 1:num_genes, FUN = function(X){
    
    gene_weights <- weights[X, ]
    
    gw <- gene_weights
    min_pos <- min(gw[gw > 0])
    if (min(gw) == 0) gw[gw == 0] <- min_pos / 2
    max_non1 <- max(gw[(1 - gw) > 1e-10])
    if (max(gw) == 1) gw[gw == 1] <- (1 + max_non1) / 2
    
    # logit transform
    response <- log(gw / (1 - gw))
    
    dat$response <- as.numeric(response)  # ensure plain numeric, no weird attributes
    
    tvals <- -1
    
    tryCatch(expr = {
      
      fit <- lme(
        fixed  = response ~ age + ethnicity + disease + ethnicity:disease + age:disease,
        random = ~ 1 | sample,
        data   = dat,
        na.action = na.omit,
        method = "REML"
      )
      
      tvals <- summary(fit)$tTable[, 4L]
      
    }, error = function(e){a <- 1}, finally = {
      
      # if error return vector of all 0's for the model coefficient t-values
      if(length(tvals) == 1)
      {
        tvals <- rep(0, 6)
        names(tvals) <- c("(Intercept)", "age", "ethnicityAfrican American",
                          "diseaseFlu-infected", "ethnicityAfrican American:diseaseFlu-infected",
                          "age:diseaseFlu-infected")
      }
      
    })
    
    return(tvals)
    
  })
  
  # MODIFY FOR YOUR VARIABLES
  int_tvals <- get_gene_t_for_term(models, "(Intercept)", genes)
  age_tvals <- get_gene_t_for_term(models, "age", genes)
  afr_tvals <- get_gene_t_for_term(models, "ethnicityAfrican American", genes)
  dis_tvals <- get_gene_t_for_term(models, "diseaseFlu-infected", genes)
  afr_dis_tvals <- get_gene_t_for_term(models, "ethnicityAfrican American:diseaseFlu-infected", genes)
  age_dis_tvals <- get_gene_t_for_term(models, "age:diseaseFlu-infected", genes)
  
  # MODIFY FOR YOUR VARIABLES
  # save to file
  save(list = c("int_tvals", "age_tvals", "afr_tvals", "dis_tvals", "afr_dis_tvals", "age_dis_tvals"), file = paste0(results_dir, '/tvals.RData'))
  
}

# run a survival analysis on mola results
# "surv_data" is a dataframe with one row per loop containing the sample-level statistics and survival data
# "weights" is a matrix with one row per loop, columns being the genes and entries being the weights (from gene_weights.RData files output from run_mola)
# "genes" is the vector of shared genes, corresponding to the rows of weights
survival_analysis <- function(surv_data, weights, genes){
  
  pvals <- c()
  for(i in 1:length(genes))
  {
    g <- genes[[i]]
    temp <- surv_data
    temp$gene_weight <- unname(weights[i, ])
    
    # MODIFY FOR YOUR ANALYSIS
    f <- as.formula(paste0(
      "Surv(PFI.time, PFI) ~ pathologic_stage + breast_carcinoma_estrogen_receptor_status + lab_proc_her2_neu_immunohistochemistry_receptor_status + gene_weight"
    ))
    
    tryCatch(expr = {
      
      # Cox regression
      res <- summary(coxph(
        formula = f,
        data = temp
      ))
      n <- nrow(res$coefficients)
      pvals <- c(pvals, res$coefficients[n, 5L])
      
    }, error = function(e){
      
      # if there was an error, add p-value 1
      pvals <- c(pvals, 1)
      
    })
    
  }
  
  signif_inds <- which(p.adjust(pvals, method = "BH") < 0.05)
  if(length(signif_inds) == 0)
  {
    return(list())
  }
  
  # refit just the models for the significant genes
  signif_genes <- genes[signif_inds]
  res <- lapply(X = signif_inds, FUN = function(X){
    
    g <- genes[[i]]
    temp <- surv_data
    temp$gene_weight <- unname(weights[i, ])
    
    # MODIFY FOR YOUR ANALYSIS
    f <- as.formula(paste0(
      "Surv(PFI.time, PFI) ~ pathologic_stage + breast_carcinoma_estrogen_receptor_status + lab_proc_her2_neu_immunohistochemistry_receptor_status + gene_weight"
    ))
    
    # Cox regression
    res <- summary(coxph(
      formula = f,
      data = temp
    ))
    
    return(res)
    
  })
  names(res) <- signif_genes
  
  return(res)
  
}



