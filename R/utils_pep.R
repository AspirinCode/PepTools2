################################################################################
# PepTools - An R-package for making immunoinformatics accessible              #
#     Copyright (C) 2017 Leon Eyrich Jessen                                    #
################################################################################
# This program is free software: you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
# This program is distributed in the hope that it will be useful,              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.       #
################################################################################




















################################################################################
################################################################################
################################################################################
#' Score a peptide against a PSSM using sum of positional scores.
#'
#' @param pep A vector of equal length peptides to be scored
#' @param PSSM A Position Specific Scoring Matrix, where the number of rows must be equal to the length of the input peptides.
#' @return A vector of scores, one for each input peptides. Scores are sorted according to input peptides.
#' @examples
#' PSSM = pssm_empty(npos = 9)
#' PSSM[1:9,1:20] = rnorm(180)
#' pep_score(c("RQGQDHPTM","RGQKTTDNA","NILYEYWDY"),PSSM)
pep_score = function(pep, PSSM){
  pep_check(pep = pep)
  pssm_check(PSSM = PSSM)
  n_pep = length(pep)
  l_pep = nchar(pep)[1]
  if( l_pep != nrow(PSSM) ){
    stop("Number of positions in peptide must match number of rows in PSSM")
  }

  # The following is tricky in order to avoid looping over millions of peptides

  # Select columns from PSSM according to all peptides concatenated
  # to one long vector, where each element is a single residue
  pep_long  = pep %>% str_c(collapse = '') %>% str_split('') %>% unlist
  pssm_long = PSSM[,pep_long]

  # Now create repeated column bound diagonal matrices and multiply with
  # long PSSM effectively setting all non-diagonal scores to zero
  # Sum of columns will now be the score of each residue at each position
  s_vec = colSums( pssm_long * matrix(rep(diag(l_pep),n_pep), nrow = l_pep) )

  # Convert to matrix per l_pep positions
  s_mat = matrix(s_vec, ncol = l_pep, byrow = TRUE)

  # and compute pep sum scores as sum of rows
  scores = rowSums(s_mat)

  # Done
  return(scores)
}
################################################################################
################################################################################
################################################################################









################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################










################################################################################
################################################################################
################################################################################
#' Translate all non-standard amino acid residue characters with 'X'
#' @param pep A character vector of peptides to be cleaned
#' @return A character vector of cleaned peptides
#' @examples
#' pep_clean(sample(letters,20,replace=TRUE))
pep_clean = function(pep){
  return( pep %>% toupper %>%
            str_replace_all(pattern = "[^ARNDCQEGHILKMFPSTWYVX-]", replacement = "X") )
}
################################################################################
################################################################################
################################################################################










################################################################################
################################################################################
################################################################################
#' Get peptide images
#'
#' Peptides are encoded using \code{pep_encode()}, resulting in peptide 'images'
#'
#' Each position in the peptide become a vector of 20 values, corresponding to
#' the rounded log odds ratio for substituting the amino acid in the peptide
#' with each of the 20 standard proteogenic amino acids. The final result is a
#' list of peptide 'images'
#'
#' @param pep A character vector of peptides to be converted to 'images'
#' @return A list of peptide 'images'
#' @examples
#' pep_get_images(pep_ran(k=9,n=10))
pep_get_images = function(pep){

  # Check input vector
  pep_check(pep = pep)

  # Set encoding matrix
  bl62_prob = BLOSUM62_PROB

  # Then we convert the vector of peptides to a matrix
  # with dimensions 'm x n' = 'n_peps x length_peps'
  p_mat = pep %>% pep_mat

  # Generate peptide 'images' and save in list
  o_list = list()
  for( i in 1:nrow(p_mat) ){
    pep_i_residues = p_mat[i,]
    pep_img = bl62_prob[pep_i_residues,]
    o_list[[i]] = pep_img
  }
  return(o_list)
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Plot peptide images
#'
#' Plot the first \code{n} encoded peptide 'images'
#'
#' Peptides are encoded using \code{pep_encode()}, resulting in peptide 'images'
#' Each position in the peptide become a vector of 20 values, corresponding to
#' the rounded log odds ratio for substituting the amino acid in the peptide
#' with each of the 20 standard proteogenic amino acids. The final result is a
#' list of peptide 'images' each with values in \code{[0;1]}.
#'
#' @param pep A character vector of peptides to be plotted as 'images'
#' @return A cowplot plot grid of the \code{n} first peptide 'images'
#' @examples
#' pep_plot_images(pep_ran(k=9,n=10))
pep_plot_images = function(pep, n = 3){

  # Check input vector
  pep_check(pep = pep)

  # Check 'n'
  if( n > length(pep) ){
    n = length(pep)
  }
  if( n < 1 ){
    stop("'n' must be larger than or equal to 1")
  }

  # Convert peptide to list of 'images'
  pep_imgs = pep[1:n] %>% pep_get_images

  # Convert each image in list to a plot
  plot_list = list()
  for( i in 1:length(pep_imgs) ){
    pep_img  = pep_imgs[[i]]
    residues = paste(rownames(pep_img), seq(1,nrow(pep_img)), sep = '_')
    residues = factor(x = residues, levels = rev(residues))
    p_img = tibble(pep_res = residues) %>%
      bind_cols(as_tibble(pep_img)) %>%
      gather(sub_res,val,-pep_res) %>%
      select(pep_res,sub_res,val) %>%
      ggplot(aes(sub_res,pep_res)) +
      geom_tile(aes(fill = val), color = "white") +
      scale_fill_gradient(low = "white", high = "black", limits = c(0,1)) +
      xlab("Score for substituting peptide residue with ...") +
      ylab("Peptide residue") +
      theme_classic() +
      theme(axis.line.x = element_blank(),
            axis.line.y = element_blank()) +
      ggtitle(paste("Peptide",pep[i],"'image' encoded"))
    plot_list[[i]] = p_img
  }

  # Done return cowplot plotlist
  return( plot_grid(plotlist = plot_list, ncol = 1) )

}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Retrieve the consensus sequence from a set of peptides
#' @param pep A vector of peptides
#' @return A string corresponding to the consensus sequence
#' @examples
#' pep_consensus(c("YMNSMQEML","FIYRHMFCV","VLFKFDMFI","KLLDRFPVA","RVLDDFTKL"))
pep_consensus = function(pep){
  get_max = function(x){
    o = names(x)[which.max(x)]
    if( length(o) > 1 ){
      o = paste(o,collapse='')
      o = paste0("[",o,"]")
    }
    return(names(x)[which.max(x)])
  }
  pep_check(pep = pep)
  pep %>% pep_mat %>% apply(2,function(p_i){ p_i %>% table %>% return }) %>%
    lapply(get_max) %>% unlist %>% paste(collapse = '') %>% return
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Encode a peptide
#'
#' Peptides are encoded using the BLOSUM62 probability matrix
#' (See BLOSUM62_PROB data matrix)
#'
#' Each position in the peptide become a vector of 20 values, corresponding to
#' the rounded log odds ratio for substituting the amino acid in the peptide
#' with each of the 20 standard proteogenic amino acids.
#'
#' The final result is a tibble with as many rows as input peptides and number
#' columns equal to the number of positions in the input peptides times the
#' number of BLOSUM62 substitions values (20) plus a column with the input
#' peptides. For 100 9-mers this yields 100 rows and 181 columns
#'
#' @param pep A character vector of peptides to be encoded
#' @return A tibble of encoded peptides
#' @examples
#' pep_encode_mat(pep_ran(k=9,n=100))
#' dim(pep_encode_mat(pep_ran(k=9,n=100)))
pep_encode_mat = function(pep){

  # Check input
  pep_check(pep)

  # Get number of positions in peptide
  n_pos = nchar(pep[1])

  # Set encoding matrix
  bl62_prob = BLOSUM62_PROB

  # Convert input peptides to one long character vector of single amino acids
  pep_str = pep %>% paste(collapse = '') %>% str_split('') %>% unlist

  # Generate output matrix, such that each position in the peptide is encoded as
  # a vector of 20 values corresponding to the values in the scaled BLOSUM62
  # matrix. Each position is then concatenated forming a vector of a length of
  # of number of positions in the peptides times the number of substitution
  # values in the BLOSUM62 matrix, i.e. for a 9-mer this will be 9 x 20 = 180.
  # The output matrix will then have dimensions 180 columns times number of input
  # peptides rows and then an extra column when includeing the input peptides
  # in the output. This way column names are:
  # A_p1 R_p1 N_p1 D_p1 C_p1 Q_p1 E_p1 G_p1 H_p1 I_p1 L_p1 K_p1 M_p1 F_p1 P_p1
  # S_p1 T_p1  W_p1  Y_p1  V_p1  A_p2 R_p2 N_p2 D_p2 C_p2 Q_p2 E_p2 G_p2 ...
  out_mat = bl62_prob[pep_str,] %>% t %>% as.vector %>%
    matrix(ncol=n_pos * ncol(bl62_prob), byrow = TRUE)
  colnames(out_mat) = paste(colnames(bl62_prob),
                            paste0('p',rep(1:n_pos,rep(ncol(bl62_prob),n_pos))),sep='_')
  out_mat = out_mat %>% as_tibble
  out_mat = as_tibble(pep) %>% bind_cols(out_mat) %>% rename(peptide=value)
  return(out_mat)
}
################################################################################
################################################################################
################################################################################
