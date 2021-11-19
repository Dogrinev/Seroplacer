#Functions:

#' Query and Reference Sequence Data Preparation
#'
#' Takes a pre-aligned set of query sequences from the user and formats the query and reference sequence file into
#' the appropriate formats for EPA-NG. This function executes the initial step of optimizing the order of the 7
#' 16S reference sequences into the best possible orientation matching the user's query. The best orientation for
#' each reference is calculated and then the reference sequences are re-sorted into the optimal orders. After this,
#' the reference sequences are all concatenated and returned in a format ready for EPA-NG. In addition, the user's
#' input query sequences are also concatenated and returned as a single FASTA line for input into EPA-NG. Note: It
#' is required that the input query sequences are pre-aligned to the main reference file using MAFFT. Non-aligned
#' sequences with varying sequence lengths will cause errors in downstream analysis steps.
#'
#' @param PathToAlignedQuery A file path leading to the directory containing the pre-aligned query sequences. This
#' file should contain a set of 7 ASVs in a FASTA format file which have been aligned to the reference sequence
#' database using MAFFT.
#'
#' @param Aligned_ReferenceFile (Provided in package) Select the correct 16S sequence reference file for the
#' organism being studied (Salmonella or E. coli).
#'
#' @return A list with 2 elements containing the user's formatted sequences. List element 1 contains the concatenated
#' query sequences. List element 2 contains the re-sorted and concatenated reference sequence file. Both of these
#' will be used as inputs for the EPA-NG wrapper function.
#'
#' @importFrom Biostrings readDNAStringSet
#'
#' @export
#'
#' @examples
#' Sequence_Outputs<-Data_Preparation(PathToAlignedQuery = "/path/to/your/file.fasta",Aligned_ReferenceFile = Full_16s_Data_Aligned)

Data_Preparation<-function(PathToAlignedQuery,Aligned_ReferenceFile)
{
  Aligned_fastaFile <- readDNAStringSet(PathToAlignedQuery,format = "fasta")
  Aligned_seq_name = names(Aligned_fastaFile)
  Aligned_sequence = paste(Aligned_fastaFile)
  Aligned_Query_FASTA <- data.frame(Aligned_seq_name, Aligned_sequence,stringsAsFactors = FALSE)
  QTestSequence<-c(Aligned_Query_FASTA[1,2],Aligned_Query_FASTA[2,2],Aligned_Query_FASTA[3,2],Aligned_Query_FASTA[4,2],Aligned_Query_FASTA[5,2],Aligned_Query_FASTA[6,2],Aligned_Query_FASTA[7,2])
  #The placer function will optimize order of references for concatenating in the next step
  Concat_Orders<-map(.x = 1:1043,.f = ~Super_HAM_Placer(TestSequence = QTestSequence,ReferenceSequence = Aligned_ReferenceFile[[.]]))
  #This function will concatenate the reference 16s alleles based on their order pulled from above
  Concatenated_References<-map(.x = 1:1043,.f = ~Concatenate_By_Order(Order = Concat_Orders[[.]],RefSequence = Aligned_ReferenceFile[[.]]))
  Prep_For_FASTA<-as.data.frame(unlist(Concatenated_References))
  Prep_For_FASTA2<-cbind(GFFs_Reference,Prep_For_FASTA)
  colnames(Prep_For_FASTA2)<-c("seq.name","seq.text")
  #Next we format and prepare prepare the query FASTA
  Query_Prep<-as.data.frame(paste(QTestSequence[1],QTestSequence[2],QTestSequence[3],QTestSequence[4],QTestSequence[5],QTestSequence[6],QTestSequence[7],sep = ""))
  Query_Prep2<-cbind("Query",Query_Prep)
  colnames(Query_Prep2)<-c("seq.name","seq.text")
  Sequence_Lists<-as.list(NULL)
  Sequence_Lists[[1]]<-Query_Prep2
  Sequence_Lists[[2]]<-Prep_For_FASTA2
  return(Sequence_Lists)
}

#' Wrapper Function for MAFFT
#'
#' Creates an R interface for using the command line tool MAFFT in an R environment. It is required that MAFFT is
#' initially installed on the machine. Installation instructions for MAC OS X and Windows can be located at
#' <https://mafft.cbrc.jp/alignment/software/source.html>. The installation path should be noted and may be needed
#' as a variable in the wrapper function. This wrapper is a more specific implementation of MAFFT which will take
#' the user's input query sequences and align them to a previously aligned reference sequence database.
#' The input sequences will be aligned with their lengths adjusted to match the reference multiple sequence
#' alignment (MSA). The input sequences will then be appended to the end of the MSA for use in downstream analysis.
#' This wrapper is not currently designed to include all options contained in MAFFT since our implementation requires
#' a specific set of parameters to be used.
#'
#' @param x A FASTA file containing the appropriate organism's multiple sequence alignment. These are pre-built and
#' provided with the package, and the correct one must be selected. The user's query will be aligned matching this
#' multiple sequence alignment.
#'
#' @param y A set of amplicon sequencing variants (ASVs) representing the 7 16S alleles. These are the query sequences
#' which will be aligned to the multiple sequence alignment in x.
#'
#' @param exec (Optional) The path on the machine in which the MAFFT program is installed. This path will default to
#' the default installation path if left blank (/usr/local/bin/mafft).
#'
#' @param options (Optional) Option settings which can be passed to the MAFFT command. The default setting will add the
#' option --keeplength. This option is required for this analysis, as it will ensure that the length of the sequences
#' in the multiple sequence alignment is not changed after the addition of the user's query sequences.
#'
#' @param file (Optional) A path to which the output FASTA should be saved. If left blank, the wrapper function will
#' output the resulting FASTA as an R environment object.
#'
#' @param add Attaches the --add argument to the MAFFT command, setting the program to append the user's query sequences
#' to the provided reference multiple sequence alignment. This command is required for merging the query sequences into
#' the reference FASTA file.
#'
#' @return A multiple sequence alignment in FASTA format containing the reference sequences and the aligned query
#' sequences appended as the last sequences in the file. The sequences in this alignment will be identical in length
#' to the original reference multiple sequence alignment.
#'
#' @importFrom ips write.fas
#' @importFrom ape read.dna
#'
#' @export
#'
#' @examples
#' Output_FASTA<-mafft_wrap(x = Full_Alignment.fasta, y = Query.fasta, options = "--keeplength")

mafft_wrap<-function(x, y, exec, options, file, add)
{
  os <- .Platform$OS
  if (missing(exec))
    exec <- "/usr/local/bin/mafft"
  if (missing(options)){
    options <- "--keeplength"
  }
  else{
    options <- match.arg(options, c("--adjustdirection",
                                    "--adjustdirectionaccurately",
                                    "--keeplength"))
    options <- paste(options, collapse = " ")
    options <- paste(rep(" ", 2), collapse = options)
  }
  if (missing(add))
    add <- "add"
  add <- match.arg(add, c("add", "addprofile"))
  add <- paste("--", add, sep = "")
  fns <- vector(length = 3)
  for (i in seq_along(fns)) fns[i] <- tempfile(pattern = "mafft",
                                               tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])
  write.fas(x, fns[1])
  write.fas(y, fns[2])
  call.mafft <- paste(exec, add, fns[2], options, fns[1],
                      ">", fns[3])
  if (os == "unix") {
    system(call.mafft, intern = FALSE, ignore.stdout = FALSE)
    res <- (file.info(fns[3])$size > 1)
    if (res != 0) {
      res <- read.dna(fns[3], format = "fasta")
    }
  }
  else {
    res <- system(call.mafft, intern = TRUE, ignore.stderr = FALSE)
    if (length(grep("error|ERROR", res))) {
      res <- 0
    }
    else {
      res <- read.dna(fns[3], format = "fasta")
    }
  }
  unlink(fns[file.exists(fns)])
  if (!missing(file)) {
    write.fas(res, file)
  }
  else {
    return(res)
  }
}

#' Wrapper function for EPA-NG
#'
#' Creates an R interface for using the command line tool EPA-NG in an R environment. It is required that EPA-NG
#' is initially installed on the machine. Installation instructions for all platforms can be found at
#' <https://github.com/Pbdas/epa-ng>. The installation path should be noted and may be needed as a variable in the
#' wrapper function. This wrapper is an implementation of EPA-NG which will take the aligned and concatenated query
#' sequences along with the aligned and concatenated reference sequence database and generate phylogenetic placement
#' results. This wrapper is not currently designed to include all options available in EPA-NG since our implementation
#' uses a limited range of required parameters.
#'
#' @param exec (Optional) The path on the machine in which the EPA-NG program is installed. This path will default to
#' the default installation path if left blank (/usr/local/bin/epa-ng).
#'
#' @param msa A reference multiple sequence alignment in FASTA format. For serovar prediction usage, this MSA should
#' be in the form of concatenated sequences.
#'
#' @param tree A phylogenetic tree representing the reference sequences in Newick format. This tree is provided as
#' part of the Seroplacer package.
#'
#' @param query A query sequence in FASTA format. For serovar prediction usage, this query sequence should be concatenated
#' and pre-aligned to the MSA.
#'
#' @param model File containing the GTRGAMMA model parameters. This model file is provided as part of the Seroplacer
#' package.
#'
#' @param filter_max (Optional) A number telling the EPA-NG function the maximum number of nodes to report in the
#' placement results. We generally want to see all of the possible placement results, if left blank this option defaults
#' to 100.
#'
#' @param file (Optional) A path to which the output jplace file should be saved. If left blank, the wrapper function
#' will output the resulting jplace as an R environment object.
#'
#' @return A JSON placement (jplace) file containing the placement results. This file will report nodes where placements
#' were made along with other relevant data, including pendant lengths and likelihood weight ratios for each placement.
#'
#' @importFrom ips write.fas
#' @importFrom ape write.tree
#' @importFrom treeio read.jplace
#'
#' @export
#'
#' @examples
#' Output.jplace<-epa_ng_wrap(msa = concatenated_references.fasta, tree = full.tree, query = query.fasta, model = info.raxml.bestModel, filter_max = 100)

epa_ng_wrap<-function(exec, msa, tree, query, model, filter_max, file)
{
  os <- .Platform$OS
  if (missing(exec))
    exec <- "/usr/local/bin/epa-ng"
  if (missing(filter_max)){
    filter_max <- c("--filter-max 100")
  }
  else{
    filter_max_option <- paste("--filter-max ",filter_max,sep="")
  }
  fns <- vector(length = 5)
  for (i in seq_along(fns)) fns[i] <- tempfile(pattern = "epang",
                                               tmpdir = tempdir(), fileext = ".txt")
  fns[5]<-tempdir()
  unlink(fns[file.exists(fns)])
  write.fas(msa,fns[1])
  write.tree(tree,fns[2])
  write.fas(query,fns[3])
  writeChar(model,fns[4],eos=NULL)
  JPlaceFile<-paste(fns[5],"/epa_result.jplace",sep="")
  call.epa_ng <- paste(exec, "--ref-msa", fns[1], "--tree",  fns[2], "--query", fns[3], "--model", fns[4], filter_max, "--redo --outdir", fns[5])
  if (os == "unix") {
    system(call.epa_ng, intern = FALSE, ignore.stdout = FALSE)
    res <- (file.info(JPlaceFile)$size > 1)
    if (res != 0) {
      res <- read.jplace(JPlaceFile)
    }
  }
  else {
    res <- system(call.epa_ng, intern = TRUE, ignore.stderr = FALSE)
    if (length(grep("error|ERROR", res))) {
      res <- 0
    }
    else {
      res <- read.jplace(JPlaceFile)
    }
  }
  unlink(fns[file.exists(fns)])
  if (!missing(file)) {
    res <- readChar(JPlaceFile,nchars = 10000000)
    writeChar(res,file,eos=NULL)
  }
  else {
    return(res)
  }
}

#' Reference Sequence Re-ordering
#'
#' Because we are working with a multi-copy marker gene, it is important to re-sort the reference 16S sequences into
#' the order which best matches the user submitted query sequences. This function will calculate the number of nucleotide
#' differences between all query and reference sequences and then determine the combination which results in the
#' minimum number of mismatches summed across all allele pairings.
#'
#' @param TestSequence A vector of nucleotide sequences containing 7 sequences representing the 7 16S query alleles.
#'
#' @param ReferenceSequence A vector of nucleotide sequences containing 7 sequences representing the 7 16s reference
#' alleles. This vector should contain 8 total elements, the first being a name, followed by 7 sequence elements.
#'
#' @return A vector of 7 numbers (1-7) which represents the order of reference sequences which contains the lowest
#' total number of nucleotide mismatches when compared to the query sequence.
#'
#' @importFrom purrr map_int
#'
#' @export
#'
#' @examples
#' Best_Sequence_Order <- Super_HAM_Plaer(TestSequence = Sequence1, ReferenceSequence = Sequence2)

Super_HAM_Placer<-function(TestSequence,ReferenceSequence)
{
  Row1<-map_int(.x = 1:7,.f = ~string.diff(a = TestSequence[.],b = ReferenceSequence[2],exclude = c("n","N","?"),ignore.case = TRUE))
  Row2<-map_int(.x = 1:7,.f = ~string.diff(a = TestSequence[.],b = ReferenceSequence[3],exclude = c("n","N","?"),ignore.case = TRUE))
  Row3<-map_int(.x = 1:7,.f = ~string.diff(a = TestSequence[.],b = ReferenceSequence[4],exclude = c("n","N","?"),ignore.case = TRUE))
  Row4<-map_int(.x = 1:7,.f = ~string.diff(a = TestSequence[.],b = ReferenceSequence[5],exclude = c("n","N","?"),ignore.case = TRUE))
  Row5<-map_int(.x = 1:7,.f = ~string.diff(a = TestSequence[.],b = ReferenceSequence[6],exclude = c("n","N","?"),ignore.case = TRUE))
  Row6<-map_int(.x = 1:7,.f = ~string.diff(a = TestSequence[.],b = ReferenceSequence[7],exclude = c("n","N","?"),ignore.case = TRUE))
  Row7<-map_int(.x = 1:7,.f = ~string.diff(a = TestSequence[.],b = ReferenceSequence[8],exclude = c("n","N","?"),ignore.case = TRUE))
  Hamming_Table_Result<-rbind(Row1,Row2,Row3,Row4,Row5,Row6,Row7)
  Combination_Integers<-map_int(.x = 1:5040,.f = ~Super_HAM_Combinator(Hamming_Table = Hamming_Table_Result,Combination = as.integer(CombinationTable[.,])))
  Best_Order<-as.integer(CombinationTable[which(Combination_Integers==min(Combination_Integers))[1],])
  return(Best_Order)
}

#' String Nucleotide Mismatch Calculation
#'
#' This function calculates the number of nucleotide position differences between two equal length character strings.
#' We use this in order to calculate all possible mismatch values between query and references sequences in order
#' to re-sort the reference sequences into the appropriate order.
#'
#' @param a A character string, equal in length to b.
#'
#' @param b A character string, equal in length to a.
#'
#' @param exclude Characters which to exclude from mismatch calculation.
#'
#' @param ignore.case TRUE/FALSE whether to ignore case of characters.
#'
#' @return A number value which is the number of nucleotide differences between the two inputted strings.
#'
#' @export
#'
#' @examples
#' string.diff(a = "AAAAACA",b = "AAAAGCA",exclude = c("n","N","?"),ignore.case = TRUE)

string.diff<-function(a,b,exclude=c("n","N","?"),ignore.case=TRUE)
{
  a<-toupper(a)
  b<-toupper(b)
  diff.a<-unlist(strsplit(a,split=""))
  diff.b<-unlist(strsplit(b,split=""))
  diff.d<-rbind(diff.a,diff.b)
  for(ex.loop in 1:length(exclude))
  {
    diff.d<-diff.d[,!(diff.d[1,]==exclude[ex.loop]|diff.d[2,]==exclude[ex.loop])]
  }
  differences<-sum(diff.d[1,]!=diff.d[2,])
  return(differences)
}

#' Combination Testing
#'
#' A function for calculating the number of mismatches in each possible combination of allele pairings between 7 query
#' alleles and 7 reference alleles. A 7x7 table is calculated with the number of nucleotide mismatches between each
#' query and reference initially. The combination function then calculates the sum of nucleotide mismatches for each
#' possible set of pairings of 7 alleles. There are a total of 7! = 5040 combinations which are tested to determine
#' which sum is the lowest, representing the best possible orientation of reference alleles relative to the query alleles.
#'
#' @param Hamming_Table A 7x7 table containing the number of nucleotide mismatches between 7 query and 7 reference alleles.
#'
#' @param Combination A vector of the numbers 1-7 in a specific combination order.
#'
#' @return An integer value which is the  sum of nucleotide mismatches for the specific 1-7 combination which was
#' tested. This function is applied over a table of all possible 5040 combinations to test all possible allele orders.
#'
#' @export
#'
#' @examples
#' Missing

Super_HAM_Combinator<-function(Hamming_Table,Combination)
{
  Sum_Output<-sum(Hamming_Table[1,Combination[1]],Hamming_Table[2,Combination[2]],Hamming_Table[3,Combination[3]],Hamming_Table[4,Combination[4]],Hamming_Table[5,Combination[5]],Hamming_Table[6,Combination[6]],Hamming_Table[7,Combination[7]])
  return(Sum_Output)
}

#' Reference Sequence Concatenation
#'
#' EPA-NG placement input requires a single sequence, so we concatenate the multi-copy alleles to resolve this. This
#' function concatenates a set of reference alleles based on the previously calculated best order determined using
#' Super_HAM_Placer.
#'
#' @param Order A vector containing the integers 1-7 in an order which the reference alleles should be sorted into for
#' matching with the query sequence.
#'
#' @param RefSequence A vector of nucleotide sequences containing 7 sequences representing the 7 16s reference
#' alleles. This vector should contain 8 total elements, the first being a name, followed by 7 sequence elements.
#'
#' @return A character string containing 7 concatenated 16S sequences, in the order provided by the order parameter.
#'
#' @export
#'
#' @examples
#' Concatenated_Sequences <- Concatenate_By_Order(Order = c(1,4,2,3,6,7,5), RefSequence = Reference_Sequence)


Concatenate_By_Order<-function(Order,RefSequence)
{
  Sequence<-RefSequence[2:8]
  Concatenated_Output<-paste(Sequence[Order[1]],Sequence[Order[2]],Sequence[Order[3]],Sequence[Order[4]],Sequence[Order[5]],Sequence[Order[6]],Sequence[Order[7]],sep = "")
  return(Concatenated_Output)
}

#' MRCA Calculation From Placement Data
#'
#' This function is the final optimized version of the clade finding algorithm which we apply to a set of placement
#' results. The steps for this algorithm are described below:
#' 1. EPA-NG is ran given the query sequence and all reported edges within the 99.9% likelihood weight ratio (LWR)
#' range are considered for possible initial hits. Query sequences are inserted in between two nodes by EPA-NG, and
#' either the distal or proximal node is chosen depending on which is closer to the query insertion point.
#' 2. All tips which are descendants of the closest node to the query branch insertion point are considered to be
#' the set of initial hits.
#' 3. The most recent common ancestor (MRCA) of these initial hits is calculated and considered the initial clade MRCA.
#' 4. The maximum pendant length from the reported edges in EPA-NG is recorded and multiplied by a scaler in order
#' to perform the pendant length adjustment. The scaler is varied in order to optimize for best algorithm performance.
#' 5. Pendant length adjustment is done by moving a pairwise distance from the initial clade MRCA towards the tree
#' root by a multiple of the maximum pendant length value. After traveling up the tree to a new position, the closest
#' node (distal or proximal) to this position is recorded and labeled as the pendant adjusted MRCA. If the pendant
#' length is small enough, it is possible that the MRCA will not change during pendant length adjustment.
#' 6. All descendants of this pendant adjusted MRCA are considered the final results of the placement.
#'
#' @param Pendant_Multi An integer which is the multiplier by which to scale the pendant length adjustment in step 4.
#'
#' @param Tree A phylogenetic tree representing the reference sequences in Newick format. This tree is provided as
#' part of the Seroplacer package.
#'
#' @return An integer which is the node value in Tree of the final pendant-adjusted MRCA.
#'
#' @importFrom treeio read.jplace
#' @importFrom ape dist.nodes
#' @importFrom purrr map
#' @importFrom phangorn Descendants
#' @importFrom tidytree rootnode
#'
#' @export
#'
#' @examples
#' MRCA <- Clade_Hit_Finder_Pendant_Final(Pendant_Multi = 1.5, Tree = full.tree)

Clade_Hit_Finder_Pendant_Final<-function(Pendant_Multi,Tree)
{
  JPlace<-read.jplace("epa_result.jplace")
  vert.tree<-Tree
  PlacedEdges<-JPlace@placements$node
  Node_Distances<-dist.nodes(x = vert.tree)
  PlacedEdges_List<-map(.x = 1:length(PlacedEdges),.f = ~Placement_Node_Selector(Distal_Node = PlacedEdges[.],Distal_Length = JPlace@placements$distal_length[.],Tree = vert.tree,Root = rootnode(vert.tree),Dist_Table = Node_Distances))
  PlacedEdges_Corrected<-unlist(PlacedEdges_List)
  All_Descendants<-map(.x = PlacedEdges_Corrected,.f = ~Descendants(x = vert.tree,node = .))
  All_Descendants_Vec<-unlist(All_Descendants)
  TipNumbersOI_Filtered<-unique(All_Descendants_Vec)
  Edge.Table<-cbind(vert.tree$edge,vert.tree$edge.length)
  Ancestor_List<-map(.x = TipNumbersOI_Filtered,.f = ~unlist(Ancestors(vert.tree,node = .)))
  Sum_List<-as.list(NULL)
  for(j in 1:length(Ancestor_List[[1]]))
  {
    Test_Vec<-map_int(.x = Ancestor_List,.f = ~match(Ancestor_List[[1]][j],.))
    Sum_List[[j]]<-Test_Vec
  }
  Result_Vec<-map_int(.x = Sum_List,.f = ~sum(.))
  Result_Vec_NoNA<-Result_Vec[!is.na(Result_Vec)]
  Target_Value<-min(Result_Vec_NoNA)
  Closest_Parent_Node<-Ancestor_List[[1]][which(Result_Vec==Target_Value)]
  Root_Node<-rootnode(vert.tree)
  if(Closest_Parent_Node == Root_Node)
  {
    return(Closest_Parent_Node)
  }
  if(Closest_Parent_Node != Root_Node)
  {
    Potential_Ancestors<-Ancestors(x = vert.tree,node = Closest_Parent_Node)
    Potential_Ancestors_Final<-c(Closest_Parent_Node,Potential_Ancestors)
    Edge_Scores<-Edge.Table[which(Edge.Table[,2] %in% Potential_Ancestors_Final),]
    Pendant<-max(JPlace@placements$pendant_length)
    Branch_Sums<-map(.x = 1:length(Edge_Scores[,3]),.f = ~sum(Edge_Scores[.:length(Edge_Scores[,3]),3]))
    Branch_Sums_Vec<-unlist(Branch_Sums)
    MRCA<-Edge_Scores[which(abs(Branch_Sums_Vec-Pendant*Pendant_Multi)==min(abs(Branch_Sums_Vec-Pendant*Pendant_Multi))),2]
    return(MRCA)
  }
}

#' Distal / Proximal Node Selection
#'
#' This is a function built into the MRCA calculation function which performs the initial determination of whether
#' each placement insertion is closer to the distal or proximal node (Step 1. from Clade_Hit_Finder_Pendant_Final).
#'
#' @param Distal_Node An integer value for the distal node resulting from a placement.
#'
#' @param Distal_Length An integer value for the distal length between the node inserted by the placement algorithm
#' and the distal node.
#'
#' @param Tree A phylogenetic tree representing the reference sequences in Newick format. This tree is provided as
#' part of the Seroplacer package.
#'
#' @param Root An integer value for the root node of the phylogenetic tree.
#'
#' @param Dist_Table A table of pairwise distances between all pairs of nodes in the phylogenetic tree. Calculated from
#' dist.nodes {ape}.
#'
#' @return An integer value representing either the distal or proximal node adjacent to the node inserted by the
#' placement algorthim
#'
#' @importFrom tidytree parent
#'
#' @export
#'
#' @examples
#' Missing

Placement_Node_Selector<-function(Distal_Node,Distal_Length,Tree,Root,Dist_Table)
{
  if(Distal_Node == Root)
  {
    return(Distal_Node)
  }
  else
  {
    Proximal_Node<-parent(.data = Tree,.node = Distal_Node)
    Total_Distance<-Dist_Table[Distal_Node,Proximal_Node]
    Proximal_Length<-Total_Distance-Distal_Length
    if(Distal_Length >= Proximal_Length)
    {
      return(Proximal_Node)
    }
    else
    {
      return(Distal_Node)
    }
  }
}

#' Serovar Prediction Results
#' This function generates the output of placement results from the final resulting clade calculated using the
#' MRCA calculation function (Clade_Hit_Finder_Pendant_Final). The MRCA calculation function determines the final
#' clade MRCA, and this function will analyze the resulting hits (all descendants of that MRCA) and report the
#' resulting serovars in a data table. This version of the results generation will trim and remove any serovars
#' below 5% representation in the final clade. For a complete view of the results, use Placement_Results_Output_Full
#' instead.
#'
#' @param MRCA An integer value for the node which is the final pendant adjusted MRCA.
#'
#' @param Tree A phylogenetic tree representing the reference sequences in Newick format. This tree is provided as
#' part of the Seroplacer package.
#'
#' @return #A data table containing the following information for each resulting serovar:
#' 1. The number of matches in the final resulting clade. 2. The fraction of total matches which that serovar
#' represents. 3. The maximum depth, or the distance from the MRCA to the farthest hit in the clade (calculated
#' for each serovar). 4. The sum of branch lengths, or the total sum of branch lengths from each hit to the MRCA
#' (calculated for each serovar).
#'
#' @importFrom phangorn Descendants
#' @importFrom purrr map
#' @importFrom ape dist.nodes
#'
#' @export
#'
#' @examples
#' Serovar_Results_Table <- Placement_Results_Output(MRCA = 1352, Tree = full.tree)

Placement_Results_Output<-function(MRCA,Tree)
{
  Descendant_List<-Descendants(x = Tree,node = MRCA)
  Descendants_Vec<-unlist(Descendant_List)
  Descendant_Assemblies<-Tree$tip.label[Descendants_Vec]
  Descendant_Serovars<-GTD_Sero_Predict_Clean2[which(GTD_Sero_Predict_Clean2$Assembly %in% Descendant_Assemblies),3]
  Descendant_Serovars_NoNA<-Descendant_Serovars[which(Descendant_Serovars!="")]
  Serovar_Report<-names(table(Descendant_Serovars_NoNA))
  Sero_Percentages<-unlist(map(.x = Serovar_Report,.f = ~length(which(Descendant_Serovars %in% .))/length(Descendant_Serovars)))
  Sero_Numbers<-unlist(map(.x = Serovar_Report,.f = ~length(which(Descendant_Serovars %in% .))))
  Node_Distances<-dist.nodes(vert.tree.correct)
  Distances<-Node_Distances[MRCA,Descendants_Vec]
  Depth_Results<-map(.x = Serovar_Report,.f = ~Depth_Calculator(Serovar = .,Dists = Distances,Serovar_Names = Descendant_Serovars))
  Depth_Results_Cols<-data.frame(t(matrix(unlist(Depth_Results),nrow=2)))
  Sero_Table<-as.data.frame(cbind(Serovar_Report,Sero_Percentages,Sero_Numbers,Depth_Results_Cols))
  colnames(Sero_Table)<-c("Serovar","Fraction","Matches in Final Clade","Maximum Depth","Sum of Branch Lengths")
  Sero_Table_Sorted <- Sero_Table[order(Sero_Table$Fraction,decreasing = TRUE),]
  Sero_Table_Cut <- Sero_Table_Sorted[which(Sero_Table_Sorted$Fraction>0.05),]
  Sero_Table_Cut$Fraction<-as.numeric(Sero_Table_Cut$Fraction)
  Sero_Table_Cut$Fraction<-Sero_Table_Cut$Fraction*100
  return(Sero_Table_Cut)
}

#' Maximum Depth and Sum of Branch Length Calculation
#'
#' This function is built into the Placement_Results_Output function which generates serovar placement results data.
#' This function calculates the maximum depth, or the distance from the MRCA to the farthest hit in the clade, and also
#' the sum of branch lengths, or the total sum of branch lengths from each hit to the MRCA.
#'
#' @param Serovar A character string containing a serovar name.
#'
#' @param Dists A subset of the table of pairwise distances between all pairs of nodes in the phylogenetic tree.
#' Calculated from dist.nodes {ape}. This subset of pairwise distances should contain specifically the distances
#' between the MRCA and its descendants.
#'
#' @param Serovar_Names A character vector of all serovar names in the set of resulting hits.
#'
#' @return A numeric vector containing the maximum depth and sum of branch lengths measured in pairwise distance.
#'
#' @export
#'
#' @examples
#' Missing

Depth_Calculator<-function(Serovar,Dists,Serovar_Names)
{
  Cut_Dists<-as.numeric(Dists[which(Serovar_Names %in% Serovar)])
  Depth<-max(Cut_Dists)
  Sum_Branches<-sum(Cut_Dists)
  Output<-c(Depth,Sum_Branches)
  return(Output)
}

#' Serovar Prediction Results
#' This function generates the output of placement results from the final resulting clade calculated using the
#' MRCA calculation function (Clade_Hit_Finder_Pendant_Final). The MRCA calculation function determines the final
#' clade MRCA, and this function will analyze the resulting hits (all descendants of that MRCA) and report the
#' resulting serovars in a data table. This version of the results generation will not trim any serovars with low
#' representation. This will report all serovars for any hits found below the final pendant adjusted MRCA.
#'
#' @param MRCA An integer value for the node which is the final pendant adjusted MRCA.
#'
#' @param Tree A phylogenetic tree representing the reference sequences in Newick format. This tree is provided as
#' part of the Seroplacer package.
#'
#' @return #A data table containing the following information for each resulting serovar:
#' 1. The number of matches in the final resulting clade. 2. The fraction of total matches which that serovar
#' represents. 3. The maximum depth, or the distance from the MRCA to the farthest hit in the clade (calculated
#' for each serovar). 4. The sum of branch lengths, or the total sum of branch lengths from each hit to the MRCA
#' (calculated for each serovar).
#'
#' @importFrom phangorn Descendants
#' @importFrom purrr map
#' @importFrom ape dist.nodes
#'
#' @export
#'
#' @examples
#' Serovar_Results_Table_Full <- Placement_Results_Output(MRCA = 1352, Tree = full.tree)

Placement_Results_Output_Full<-function(MRCA,Tree)
{
  Descendant_List<-Descendants(x = Tree,node = MRCA)
  Descendants_Vec<-unlist(Descendant_List)
  Descendant_Assemblies<-Tree$tip.label[Descendants_Vec]
  Descendant_Serovars<-GTD_Sero_Predict_Clean2[which(GTD_Sero_Predict_Clean2$Assembly %in% Descendant_Assemblies),3]
  Descendant_Serovars_NoNA<-Descendant_Serovars[which(Descendant_Serovars!="")]
  Serovar_Report<-names(table(Descendant_Serovars_NoNA))
  Sero_Percentages<-unlist(map(.x = Serovar_Report,.f = ~length(which(Descendant_Serovars %in% .))/length(Descendant_Serovars)))
  Sero_Numbers<-unlist(map(.x = Serovar_Report,.f = ~length(which(Descendant_Serovars %in% .))))
  Node_Distances<-dist.nodes(vert.tree.correct)
  Distances<-Node_Distances[MRCA,Descendants_Vec]
  Depth_Results<-map(.x = Serovar_Report,.f = ~Depth_Calculator(Serovar = .,Dists = Distances,Serovar_Names = Descendant_Serovars))
  Depth_Results_Cols<-data.frame(t(matrix(unlist(Depth_Results),nrow=2)))
  Sero_Table<-as.data.frame(cbind(Serovar_Report,Sero_Percentages,Sero_Numbers,Depth_Results_Cols))
  colnames(Sero_Table)<-c("Serovar","Fraction","Matches in Final Clade","Maximum Depth","Sum of Branch Lengths")
  Sero_Table_Sorted <- Sero_Table[order(Sero_Table$Fraction,decreasing = TRUE),]
  Sero_Table_Sorted$Fraction<-as.numeric(Sero_Table_Sorted$Fraction)
  Sero_Table_Sorted$Fraction<-Sero_Table_Sorted$Fraction*100
  return(Sero_Table_Sorted)
}

#' Simple Serovar Prediction
#'
#' This function will report a predicted serovar from Placement_Results_Output. This function requires that a serovar shows at least
#' 30% representation to report, and will report all serovars over the 30% representation threshold. If there are no serovars with
#' at least 30% presence then it will report that no serovar match was found. It is highly reccomended to analyze the results table
#' in depth to understand the accuracy of the serovar prediction, as the table will contain a complete visualization of all results.
#'
#' @param Sero A serovar results table which is the output of Placement_Results_Output.
#'
#' @return A serovar prediction if a serovar represents at least 30% of the resulting hits, otherwise no result.
#'
#' @export
#'
#' @examples
#' Missing

Sero_Result<-function(Sero)
{
  if(Sero$Fraction[1]<30)
  {
    return(paste("No Serovar Match Found"))
  }
  if(Sero$Fraction[1]>30)
  {
    Sums<-unlist(map(.x = 1:length(Sero$Fraction),.f = ~sum(Sero$Fraction[1:.])))
    Sum_Index<-which(Sums>70)[1]
    Names<-c(Sero$Serovar[1:Sum_Index])
    return(paste(Names))
  }
}

#' Legend Removal for Plotting
#'
#' This function removes the legend from the phylogeny plot for visualization.
#'
#' @return N/A
#'
#' @importFrom ggplot2 theme
#'
#' @export
#'
#' @examples
#' None

no_legend <- function() {theme(legend.position="none")}

#' Phylogeny Plotting
#'
#' This function will take the final calculated clade MRCA and generate a circular phylogeny with red lines representing
#' the path from the final clade MRCA to all of its descendants (the final hits). Because previous calculations are run
#' using an unrooted tree, this function recalculates the MRCA in the rooted tree and then colors and plots based on
#' the structure of the rooted tree for improved visualization. Note that this function will not work correctly on
#' the unrooted tree and it is important to not confuse the two different trees becasue of small node differences
#' between them. After running the phylogeny plotting function, plot the outputted object to generate the actual plot.
#' Because of the size of the actual phylogeny, it is reccomended  to save the resulting plot as a PNG in dimensions
#' of at least 15,000 x 15,000 pixels. An outside viewer should be used to open the PNG and zoom in to analyze the
#' resulting phylogeny figure. Larger dimensions will increase clarity and allow better zooming at the cost of larger
#' image size.
#'
#' @param MRCA An integer value for the node which is the final pendant adjusted MRCA.
#'
#' @param InitialTable A table used for coloring branches, provided as part of the Seroplacer package.
#'
#' @param Tree A phylogenetic tree representing the reference sequences in Newick format. This tree is provided as
#' part of the Seroplacer package. This tree should be the original UNROOTED tree.
#'
#' @return A phylogeny plot with the MRCA and descendants colored in red, representing the final clade of resulting hits.
#'
#' @importFrom purrr map
#' @importFrom phangorn Ancestors
#' @importFrom phangorn Descendants
#' @importFrom ggtree ggtree
#' @importFrom ggtree %<+%
#' @importFrom ggtree geom_tiplab
#'
#' @export
#'
#' @examples
#' Plot <- Phylogeny_Plotting(MRCA = 1352, InitialTable = ColoringTable, Tree = full.tree)

Phylogeny_Plotting<-function(MRCA,InitialTable,Tree)
{
  Hits<-unlist(Descendants(x = Tree,node = MRCA))
  TipNames<-vert.tree.correct$tip.label[Hits]
  Rooted_TipNumbers<-which(rooted.vert.tree$tip.label %in% TipNames)
  Ancestor_List<-map(.x = Rooted_TipNumbers,.f = ~unlist(Ancestors(rooted.vert.tree,node = .)))
  Sum_List<-as.list(NULL)
  for(j in 1:length(Ancestor_List[[1]]))
  {
    Test_Vec<-map_int(.x = Ancestor_List,.f = ~match(Ancestor_List[[1]][j],.))
    Sum_List[[j]]<-Test_Vec
  }
  Result_Vec<-map_int(.x = Sum_List,.f = ~sum(.))
  Result_Vec_NoNA<-Result_Vec[!is.na(Result_Vec)]
  Target_Value<-min(Result_Vec_NoNA)
  Closest_Parent_Node<-Ancestor_List[[1]][which(Result_Vec==Target_Value)]
  EdgesToColor<-Descendants(x=rooted.vert.tree,node = Closest_Parent_Node,type = "all")
  InitialTable[which(InitialTable[,1] %in% EdgesToColor),2]<-"red"
  InitialTable$node<-as.integer(InitialTable$node)
  p <- ggtree(rooted.vert.tree,layout = "circular") +
    xlim(-0.5,5) +
    no_legend()
  p2 <- p %<+% TaxID_Table3 +
    geom_tiplab(aes(fill = factor(Serovar)),
                color = "black", # color for label font
                geom = "label",  # labels not text
                label.padding = unit(0.10, "lines"), # amount of padding around the labels
                label.size = 0.01,
                size = 2)
  p3 <- p2 %<+% InitialTable + aes(color=I(color))
  return(p3)
}
