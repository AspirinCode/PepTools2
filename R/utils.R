# ------------------------------------------------------------------------------
# This R-script file contain all the PepTools support functions
# ------------------------------------------------------------------------------

# Easy harvest of HTML tables
# Thu Aug  2 11:15:53 2018 ------------------------------
.get_html_table = function(x, table_no = NULL){
  library('rvest')
  my_url   = x
  my_html  = read_html(x = my_url)
  my_table = html_table(x = my_html, fill = TRUE)
  if( !is.null(table_no) ){ my_table = my_table[[table_no]] }
  return(my_table)
}

# Set background frequencies
# Thu Aug  2 11:14:52 2018 ------------------------------
# Run as .set_bg_freqs()
.set_bg_freqs = function(){
  source_url = "http://isoelectricpointdb.org/statistics.html"
  BGFREQS = .get_html_table(x = source_url, table_no = 3)
  rownames(BGFREQS) = BGFREQS$Kingdom
  BGFREQS = BGFREQS[,AMINOACIDS$three]
  aa_one = aa_translate(x = colnames(BGFREQS), to = 'one')
  colnames(BGFREQS) = aa_one
  BGFREQS = BGFREQS / 100
  save(BGFREQS, file = "data/BGFREQS.RData")
  return(0)
}

# Set a data tibble with names and abbreviations of the standard 20
# proteogenic amino acid residues
# Thu Aug  2 11:16:38 2018 ------------------------------
# Run as .set_aminoacids()
.set_aminoacids = function(){
  AMINOACIDS = data.frame(
    full      = c('Alanine', 'Arginine', 'Asparagine', 'Aspartate', 'Cysteine',
                  'Glutamine', 'Glutamate', 'Glycine', 'Histidine',
                  'Isoleucine', 'Leucine', 'Lysine', 'Methionine',
                  'Phenylalanine', 'Proline', 'Serine', 'Threonine',
                  'Tryptophan', 'Tyrosine', 'Valine'),
    three     = c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His',
                  'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp',
                  'Tyr', 'Val'),
    one       = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                  'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
    chemistry = c('Hydrophobic', 'Basic', 'Neutral', 'Acidic', 'Polar',
                  'Neutral', 'Acidic', 'Polar', 'Basic', 'Hydrophobic',
                  'Hydrophobic', 'Basic', 'Hydrophobic', 'Hydrophobic',
                  'Hydrophobic', 'Polar', 'Polar', 'Hydrophobic', 'Polar',
                  'Hydrophobic'),
    stringsAsFactors = FALSE
  )
  save(AMINOACIDS, file = "data/AMINOACIDS.RData")
  return(0)
}

# Set sample peptides
# Thu Aug  2 11:17:44 2018 ------------------------------
# Run as .set_peptides()
.set_peptides = function(){
  set.seed(509279)
  p_file = paste0("https://raw.githubusercontent.com/leonjessen/",
                  "keras_tensorflow_demo/master/data/",
                  "ran_peps_netMHCpan40_predicted_A0201_reduced_cleaned_balanced.tsv")
  PEPTIDES = read.table(file = p_file, header = TRUE, stringsAsFactors = FALSE)
  PEPTIDES = PEPTIDES[PEPTIDES$label_chr=="SB",]
  PEPTIDES = PEPTIDES[sample(1:5000),]
  PEPTIDES = PEPTIDES$peptide
  save(PEPTIDES, file = "data/PEPTIDES.RData")
  return(0)
}

# BLOSUM matrices are available via the NCBI ftp site
# Thu Aug  2 11:21:57 2018 ------------------------------
# Run as .get_BLOSUM62_matrix()
.get_BLOSUM62_matrix = function(){

  # Set ftp
  bl62_ftp = 'ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62'

  # Download as data frame
  BLOSUM62 = read.table(file = bl62_ftp, comment = '#')

  # Set and extract the 20 proteogenic amino acids
  proteogenic_aa = c('A','R','N','D','C','Q','E','G','H','I',
                     'L','K','M','F','P','S','T','W','Y','V')
  BLOSUM62 = BLOSUM62[which(rownames(BLOSUM62) %in% proteogenic_aa),
                      which(rownames(BLOSUM62) %in% proteogenic_aa)]

  # Save data
  save(BLOSUM62, file = "data/BLOSUM62.RData")

  # Set encoding matrix
  BLOSUM62_enc = BLOSUM62 / 5

  # Add null vectors for encoding X
  BLOSUM62_enc['X',] = rep(0, ncol(BLOSUM62_enc))
  BLOSUM62_enc[,'X'] = rep(0, nrow(BLOSUM62_enc))
  BLOSUM62_enc       = as.matrix(BLOSUM62_enc)

  # Save data
  save(BLOSUM62_enc, file = "data/BLOSUM62_enc.RData")

  # Done
  return(0)
}

# BLOSUM matrices are available via the NCBI ftp site
# Thu Aug  2 11:22:32 2018 ------------------------------
# Run as .get_BLOSUM50_matrix()
.get_BLOSUM50_matrix = function(){

  # Set ftp
  bl50_ftp = 'ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM50'

  # Download as data frame
  BLOSUM50 = read.table(file = bl50_ftp, comment = '#')

  # Set and extract the 20 proteogenic amino acids
  proteogenic_aa = c('A','R','N','D','C','Q','E','G','H','I',
                     'L','K','M','F','P','S','T','W','Y','V')
  BLOSUM50 = BLOSUM50[which(rownames(BLOSUM50) %in% proteogenic_aa),
                      which(rownames(BLOSUM50) %in% proteogenic_aa)]

  # Save data
  save(BLOSUM50, file = "data/BLOSUM50.RData")

  # Set encoding matrix
  BLOSUM50_enc = BLOSUM50 / 5

  # Add null vectors for encoding X
  BLOSUM50_enc['X',] = rep(0, ncol(BLOSUM50_enc))
  BLOSUM50_enc[,'X'] = rep(0, nrow(BLOSUM50_enc))
  BLOSUM50_enc       = as.matrix(BLOSUM50_enc)

  # Save data
  save(BLOSUM50_enc, file = "data/BLOSUM50_enc.RData")

  # Done
  return(0)
}

# Based on the BLOSUM50 matrix, do a PCA aiming at reducing number of
# encoding values needed per amino acid
# Tue Dec  4 11:30:32 2018 ------------------------------
# Run as .set_BLOSUM50_pca_matrix()
.set_BLOSUM50_pca_matrix = function(){

  # Add 'X' row
  BLOSUM50_X = BLOSUM50
  BLOSUM50_X['X',] = rep(0, ncol(BLOSUM50_X))

  # Do PCA and get and set rotated matrix
  BLOSUM50_pca = prcomp(x = BLOSUM50_X, center = TRUE, scale. = TRUE)$x

  # Save data
  save(BLOSUM50_pca, file = "data/BLOSUM50_pca.RData")

  # Done
  return(0)
}

# Based on the BLOSUM62 matrix, do a PCA aiming at reducing number of
# encoding values needed per amino acid
# Tue Dec  4 11:31:07 2018 ------------------------------
# Run as .set_BLOSUM62_pca_matrix()
.set_BLOSUM62_pca_matrix = function(){

  # Add 'X' row
  BLOSUM62_X = BLOSUM62
  BLOSUM62_X['X',] = rep(0, ncol(BLOSUM62_X))

  # Do PCA and get and set rotated matrix
  BLOSUM62_pca = prcomp(x = BLOSUM62_X, center = TRUE, scale. = TRUE)$x

  # Save data
  save(BLOSUM62_pca, file = "data/BLOSUM62_pca.RData")

  # Done
  return(0)
}
