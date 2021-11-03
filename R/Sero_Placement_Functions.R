#Functions:

#The data preparation function will take input data files including raw query ASVs and aligned query ASVs and generate
#properly aligned files in the working directory for further use in EPA-NG. The raw query fasta should be provided
#by the user from the data they want to place, and the aligned query fasta should be generated first using MAFFT
#before running the data preparation. The ReferenceFile and Aligned_ReferenceFile are provided and included in 
#the required files. 

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

#Wrapper function for MAFFT to enable use from r of command line program. 
mafft_wrap<-function(x, y, exec, options, file, add)
{
  os <- .Platform$OS
  if (missing(exec)) 
    exec <- "/usr/local/bin/mafft"
  if (missing(options)){
    options <- " "    
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

#Wrapper function for EPA-NG to enable use from R for this command line program - still incomplete
epa_ng_wrap<-function(exec, msa, tree, query, model, filter_max)
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
  unlink(fns[file.exists(fns)])
  write.fas(msa,fns[1])
  fns[2]<-tree
  fns[3]<-query
  fns[4]<-model
  call.epa_ng <- paste(exec, add, fns[2], options, fns[1], ">", fns[3])
  
}

#This function calculates the best orientation of ASVs to order for use in the placement tool. Since we do not know 
#which ASVs match best with which reference 16s sequences, we need to determine the best pairings. This function will
#calculate all possible combinations and determine which order is the best match, then it will output the best possible
#sequence order for each reference sequence. 

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

#This function calculates the number of nucleotide position differences between two equal length strings in order
#to determine the best order of references to align. This function is a replacement for the original use of nwhamming
#in order to work entirely with pre-aligned files. Using this function saves ~10 minutes of calculations since we 
#can avoid aligning twice. 

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

#The combinator function matches up the results in all possible combinations for finding the minimum value, this function is used
#in Super_HAM_Placer. 

Super_HAM_Combinator<-function(Hamming_Table,Combination)
{
  Sum_Output<-sum(Hamming_Table[1,Combination[1]],Hamming_Table[2,Combination[2]],Hamming_Table[3,Combination[3]],Hamming_Table[4,Combination[4]],Hamming_Table[5,Combination[5]],Hamming_Table[6,Combination[6]],Hamming_Table[7,Combination[7]])
  return(Sum_Output)
}

#This function concatenates reference sequences based on the previously determined concatenation order from Super_HAM_Placer.

Concatenate_By_Order<-function(Order,RefSequence)
{
  Sequence<-RefSequence[2:8]
  Concatenated_Output<-paste(Sequence[Order[1]],Sequence[Order[2]],Sequence[Order[3]],Sequence[Order[4]],Sequence[Order[5]],Sequence[Order[6]],Sequence[Order[7]],sep = "")
  return(Concatenated_Output)
}

#This function is the final optimized version of the clade finding algorithm which we apply to a set of placement results. 
#Description of Pendant Adjusted Placement Method v4:

#1. EPA-NG is ran given the query sequence and all reported edges within the 99.9% likelihood weight ratio (LWR) range are 
#considered for possible initial hits. Query sequences are inserted in between two nodes by EPA-NG, and either the distal 
#or proximal node is chosen depending on which is closer to the query insertion point.

#2. All tips which are descendants of the closest node to the query branch insertion point are considered to be the set of 
#initial hits.

#3. The MRCA of these initial hits is calculated and considered the initial clade MRCA.

#4. The maximum pendant length from the reported edges in EPA-NG is recorded and multiplied by a scaler in order to perform 
#the pendant length adjustment. The scaler is varied in order to optimize for best algorithm performance.

#5. Pendant length adjustment is done by moving a pairwise distance from the initial clade MRCA towards the tree root by a 
#multiple of the maximum pendant length value. After traveling up the tree to a new position, the closest node (distal or 
#proximal) to this position is recorded and labeled as the pendant adjusted MRCA. If the pendant length is small enough, it 
#is possible that the MRCA will not change during pendant length adjustment.

#6. All descendants of this pendant adjusted MRCA are considered the final results of the placement.

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

#Placement_Node_Selector is a function built into the Clade Hit Finding algorithm which performs the initial determination of 
#whether each placement is closer to a distal or proximal node (Step 1 from Clade_Hit_Finder_Pendant_Final).

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

#This function generates the output of placement results from the final resulting clade calculated using the previous function.
#The clade hit finder determines the final clade MRCA, and this function will analyze the resulting hits (all descendants of 
#that MRCA) and report the resulting serovars in a table. 
#This function will output a table containing the following information for each resulting serovar: 1. The number of matches 
#in the final resulting clade 2. The fraction of total matches which that serovar represents 3. The maximum depth,
#or the distance from the MRCA to the farthest hit in the clade (calculated for each serovar) 4. The sum of branch 
#lengths, or the total sum of branch lengths from each hit to the MRCA (calculated for each serovar)

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

#The Depth_Calculator is a function built into Placement_Results_Output for calculating the maximum depth and sum of branch
#lengths for each serovar result. 

Depth_Calculator<-function(Serovar,Dists,Serovar_Names)
{
  Cut_Dists<-as.numeric(Dists[which(Serovar_Names %in% Serovar)])
  Depth<-max(Cut_Dists)
  Sum_Branches<-sum(Cut_Dists)
  Output<-c(Depth,Sum_Branches)
  return(Output)
}

#This is a variation of Placement_Results_Output from above which will not trim the resulting data table for the user. The normal
#version of this function cuts any serovars below 5% representation for neatness in the final data reported. For users who want 
#to see the complete results including low representation serovars, this function should be used to observe the final data. 

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

#This is a simple function to report a predicted serovar. It requires that a serovar shows at least 30% representation to report,
#otherwise it will report that no serovar match was found. 

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

#This is a simple plotting function built into Phylogeny_Plotting to remove the legend for plotting.

no_legend <- function() theme(legend.position="none")

#Phylogeny_Plotting will take the final calculated clade MRCA and generate a circular phylogeny with red lines representing the 
#path from the final clade MRCA to all of its descendants (the final hits). Because previous calculations are run using an 
#unrooted tree, this function recalculates the MRCA in the rooted tree and then colors and plots based on the structure of the 
#rooted tree for improved visualization. Note that this function will not work correctly on the unrooted tree and it is important
#to not confuse the two different trees becasue of small node differences between them. 

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