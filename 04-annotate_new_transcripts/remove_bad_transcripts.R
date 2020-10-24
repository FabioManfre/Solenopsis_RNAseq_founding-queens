#Script for removing transcripts exceeding the length of the scaffold in the reference genome:
#Load the data into R:

#.gtf file
annot <- read.table("transcripts.gtf", sep = "\t")
colnames(annot) <- c("Scaffold", "Origin", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute")

#.fai file
index <- read.table("Si_gnG.fa.fai")
colnames(index) <- c("Scaffold", "Length", "Offset", "Linebases", "Linewidth")

#The 5th column of the gtf file shows the end of a feature in the .gtf file. This number should never
#be longer than that on the 2nd column of the index file: the length of the scaffold in the gnG assembly.

#I will focus on the whole transcript only, and then remove, when necessary, the exons associated to it:
transcripts_only <- annot[annot$Feature == "transcript", ]

#Merge the two files to make the scaffolds coincide with their transcript:
#Only the scaffolds which contain transcripts will be retained

merged_scaffolds <- merge(transcripts_only, index)

#Remove those transcripts which exceed the length of the scaffold
#Create a vector to store the position of those transcripts:
to_remove <- NA
for (i in 1:nrow(merged_scaffolds)){
	if(merged_scaffolds$End[i] > merged_scaffolds$Length[i]){
		to_remove[length(to_remove) + 1] <- i
	}
}
to_remove <- to_remove[-1]

#Sort the transcripts list so that the order is the same as in merged_scaffolds:
transcripts_only <- transcripts_only[order(transcripts_only$Scaffold), ]
#Sanity check

if(!identical(transcripts_only$Scaffold, merged_scaffolds$Scaffold)){ warning("merged_scaffolds and trnascripts_only do not have the same scaffolds in it")}

#Store the bad transcripts:
bad_transcripts  <- transcripts_only[to_remove, ]

#Retrieve the transcript id from the bad transcripts:

bad_transcripts$Attribute <- as.character(bad_transcripts$Attribute)
bad_transcript_id 	  <- strsplit(bad_transcripts$Attribute, ";") 
bad_transcript_id 	  <- unlist(lapply(bad_transcript_id, '[[', 2))

#Identify and remove the bad transcripts and their exons from the annotation file:
bad_features_id <- grep(paste(bad_transcript_id, collapse="|"), annot$Attribute)
annot_clean <- annot[-bad_features_id, ]

#Identifiy and remove features with empty gene id_attribute
empty_gene_id <- grep("gene_id ;", annot_clean$Attribute)
annot_clean <- annot_clean[-empty_gene_id,]

#Output the clean annotation file
write.table(annot_clean, "transcripts_clean.gtf", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)











 

