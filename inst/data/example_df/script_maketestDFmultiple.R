# make test datasets for cgmappeR

# df: multiple samples/patients per sample group
a <- epic450anno; a <- a[a$pos>=7589868 & a$pos<=7591868 & a$chr=="chr17",]
cgs <- a$Name[1:20]
samp1a <- sample(100,20)/100; samp2a <- sample(100,20)/100; samp3a <- sample(100,20)/100
samp1b <- sample(100,20)/100; samp2b <- sample(100,20)/100; samp3b <- sample(100,20)/100

dfi.3s <- data.frame(cpg=rep(cgs,6),sampGroups=c(rep("A",60),rep("B",60)),
                     methyl=c(samp1a,samp2a,samp3a,samp1b,samp2b,samp3b),
                     sampID=c(rep("1",20),rep("2",20),rep("3",20),rep("1",20),rep("2",20),rep("3",20)),
                     stringsAsFactors = F)
min(a[c(1:20),]$pos) # 7589976
max(a[c(1:20),]$pos) # 7591838

write.csv(dfi.3s,file="methyldf_test_2grp3samp.csv",row.names = F)