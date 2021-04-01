
result$level = NA
result[result$Name == "LUCA|1",]$level = 1
result[result$Name == "Eukaryota|1", ]$level = 2
result[result$Name == "Opisthokonta|1", ]$level = 3
result[result$Name == "Metazoa|1", ]$level = 4
result[result$Name == "Eumetazoa|1", ]$level = 5
result[result$Name == "Bilateria|1", ]$level = 6 
result[result$Name == "Deuterostomia|1", ]$level = 7
result[result$Name == "Chordata|1", ]$level = 8
result[result$Name == "Vertebrata|1", ]$level = 9
result[result$Name == "Euteleostomi|1", ]$level = 10
result[result$Name == "Sarcopterygii|1", ]$level = 11
result[result$Name == "Tetrapoda|1", ]$level = 12
result[result$Name == "Amniota|1", ]$level = 13
result[result$Name == "Mammalia|1", ]$level = 14 
result[result$Name == "Theria|1", ]$level = 15
result[result$Name == "Eutheria|1", ]$level = 16
result[result$Name == "Boreoeutheria|1", ]$level = 17
result[result$Name == "Euarchontoglires|1", ]$level = 18
result[result$Name == "Primates|1", ]$level = 19
result[result$Name == "Haplorrhini|1", ]$level = 20
result[result$Name == "Simiiformes|1", ]$level = 21
result[result$Name == "Catarrhini|1", ]$level = 22
result[result$Name == "Hominoidea|1", ]$level = 23
result[result$Name == "Hominidae|1",]$level = 24
result[result$Name == "Homininae|1",]$level = 25

clade_gene = table(result$level)
par(mar=c(7,6,1,2))
barplot(clade_gene,col = 'black',xlab = 'Phylostrata level',ylab = 'Number of genes')
