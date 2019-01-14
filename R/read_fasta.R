#' Read a FASTA file
#'
#' read_fasta reads a FASTA file of nucleotides or amino acids file and returns
#' a tibble with number of rows corresponding to the number of sequences and two
#' variables: 'fasta_header' and 'sequences'
#' @param file The FASTA file to be read either local or URL
#' @examples
#' read_fasta(file = 'my_fasta_file.fsa')
#' read_fasta(file = 'https://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt')
#'
#' @export
read_fasta = function(file){

  # Check arguments
  if( !is.character(file) ){
    stop("'file' has to be a string specifying a file name")
  }
  if( length(file) > 1 ){
    stop("Only one file name")
  }
  if( !file.exists(file) ){
    stop(paste("Unable to read file",file))
  }

  # Read lines from file
  file_connection = file(file, "r")
  lines = readLines(con = file_connection)
  close(file_connection)

  # Remove any comments
  keep  = which(substr(lines, 1, 1) != '#')
  lines = lines[keep]

  # Set variables
  n_seqs    = sum(substr(lines, 1, 1) == '>')
  headers   = rep(NA, n_seqs)
  sequences = rep('', n_seqs)

  # Iterate over lines and fill containers
  entry_no = 1
  for( line in lines ){
    if( grepl(pattern = "^>", x = line) ){
      headers[entry_no] = line
      entry_no = entry_no + 1
    } else {
      sequences[entry_no-1] = paste0(sequences[entry_no-1], line)
    }
  }

  # Delete lines
  rm(lines)

  # Return output tibble
  return(data.frame(fasta_header = headers, sequence = sequences,
                    stringsAsFactors = FALSE))
}
