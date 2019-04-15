library(ggseqlogo)
library(rlist)




dnvs = read.table("~/Desktop/MNV/180221_dnvs_allcombined_context.tsv",sep="\t",header=T)
nucs = c("A", "C", "G", "T")
dinuc_pairs = c()
for (i in nucs){
  for (j in nucs){
    dinuc_pairs = c(dinuc_pairs, paste0(i,",",j))
  }
}
for (refs in dinuc_pairs){
  l = list()
  lnames = c()
  for (alts in dinuc_pairs){
    if (substring(refs,1,1)!=substring(alts,1,1) & substring(refs,3,3)!=substring(alts,3,3)){    #nullじゃないのは9通り only.
      patterns = as.vector(dnvs$localseq_22b[dnvs$refs==refs & dnvs$alts==alts & (dnvs$pass_num==2)])
      l = list.append(l, patterns)
      lnames = c(lnames, paste0(refs, "to", alts, ", n=", length(patterns)))
    }
  }
  names(l) = lnames
  png(paste0(figdir, "ref", refs, "_bits_logo.png"), width = 1200, height = 400, res=108)
  print (ggseqlogo(l, method='bits'))
  dev.off()
}
for (refs in dinuc_pairs){
  l = list()
  lnames = c()
  for (alts in dinuc_pairs){
    if (substring(refs,1,1)!=substring(alts,1,1) & substring(refs,3,3)!=substring(alts,3,3)){    #nullじゃないのは9通り only.
      patterns = as.vector(dnvs$localseq_22b[dnvs$refs==refs & dnvs$alts==alts])
      l = list.append(l, patterns)
      lnames = c(lnames, paste0(refs, "to", alts, ", n=", length(patterns)))
    }
  }
  names(l) = lnames
  png(paste0(figdir, "filteredout_ref", refs, "_bits_logo.png"), width = 1200, height = 400, res=108)
  print (ggseqlogo(l, method='bits'))
  dev.off()
}


#logo specific for those repetitive ones:
#and also for the standardized ones
d1 = read.table("~/Downloads/chr1_context_std.tsv",sep="\t",header=T)
#first, simply collapsing the revcomp
d1$refs = paste0(substring(d1$ref_std,1,1), ",",substring(d1$ref_std,2,2))
d1$alts = paste0(substring(d1$alt_std,1,1), ",",substring(d1$alt_std,2,2))
nucs = c("A", "C", "G", "T")
dinuc_pairs = c()
for (i in nucs){
  for (j in nucs){
    dinuc_pairs = c(dinuc_pairs, paste0(i,",",j))
  }
}
for (refs in dinuc_pairs){
  l = list()
  lnames = c()
  for (alts in dinuc_pairs){
    if (substring(refs,1,1)!=substring(alts,1,1) & substring(refs,3,3)!=substring(alts,3,3)){    #nullじゃないのは9通り only.
      patterns = as.vector(d1$context_ref_std[d1$refs==refs & d1$alts==alts])
      l = list.append(l, patterns)
      lnames = c(lnames, paste0(refs, "to", alts, ", n=", length(patterns)))
    }
  }
  names(l) = lnames
  png(paste0("~/Downloads/d1_ref", refs, "_bits_logo_collapsed.png"), width = 1200, height = 400, res=108)
  print (ggseqlogo(l, method='bits'))
  dev.off()
}

#second, focusing on repetitive ones
d1rep = d1[d1$is_repetitive==1,]
d1rep$refs = paste0(substring(d1rep$ref_std,1,1), ",",substring(d1rep$ref_std,2,2))
d1rep$alts = paste0(substring(d1rep$alt_std,1,1), ",",substring(d1rep$alt_std,2,2))
nucs = c("A", "C", "G", "T")
dinuc_pairs = c()
for (i in nucs){
  for (j in nucs){
    dinuc_pairs = c(dinuc_pairs, paste0(i,",",j))
  }
}
for (refs in dinuc_pairs){
  l = list()
  lnames = c()
  for (alts in dinuc_pairs){
    if (substring(refs,1,1)!=substring(alts,1,1) & substring(refs,3,3)!=substring(alts,3,3)){    #nullじゃないのは9通り only.
      patterns = as.vector(d1rep$context_ref_std[d1rep$refs==refs & d1rep$alts==alts])
      l = list.append(l, patterns)
      lnames = c(lnames, paste0(refs, "to", alts, ", n=", length(patterns)))
    }
  }
  names(l) = lnames
  png(paste0("~/Downloads/d1rep_ref", refs, "_bits_logo_collapsed.png"), width = 1200, height = 400, res=108)
  print (ggseqlogo(l, method='bits'))
  dev.off()
}


#とりあえずいくつか重要なものに関してだけbitsやる
d1rep = d1[d1$is_repetitive==1,]
d1rep$refs = paste0(substring(d1rep$ref_std,1,1), ",",substring(d1rep$ref_std,2,2))
d1rep$alts = paste0(substring(d1rep$alt_std,1,1), ",",substring(d1rep$alt_std,2,2))
refsall = c("A,C","C,A","A,T","T,A","A,A","A,A")
altsall = c("C,A","A,C","T,A","A,T","T,T","C,C")
l = list()
lnames = c()
for (i in seq(length(refsall))){
  refs = refsall[i]
  alts = altsall[i]
  patterns = as.vector(d1rep$context_ref_std[d1rep$refs==refs & d1rep$alts==alts])
  l = list.append(l, patterns)
  lnames = c(lnames, paste0(refs, "to", alts, ", n=", length(patterns)))
}
names(l) = lnames
png(paste0("~/Downloads/bits_logo_fig.png"), width = 1200, height = 400, res=108)
print (ggseqlogo(l, method='bits'))
dev.off()


#ggseqlogoのperformance evaluation. どこまでいけるか.

#focus on repeats: AC->CA, CA->AC, AT->TA, TA->AT, AA->TT, AA->CC for d=1
#also annotated ns manually
ns = c(8401, 8688, 5691, 7139, 27616, 7160)
for (ptn in c("AC->CA", "CA->AC", "AT->TA", "TA->AT", "AA->TT", "AA->CC")){
  t = read.table(paste0("~/Downloads/",ptn, "_d1.tsv"), sep="\t")
  t = as.character(t$V1) #make it a vector
  p1 = ggseqlogo(t)
  print(p1)
}
#OK now plot these in a single figure
l = list()
ns = c(8401, 8688, 5691, 7139, 27616, 7160)
ptns = c("AC->CA", "CA->AC", "AT->TA", "TA->AT", "AA->TT", "AA->CC")
for (ptn in ptns){
  t = read.table(paste0("~/Downloads/",ptn, "_d1.tsv"), sep="\t")
  t = as.character(t$V1) #make it a vector
  l = list.append(l,t)
}
names(l) = paste0(ptns, ", n=", ns)
png(paste0("~/Downloads/bits_logo_fig3_final.png"), width = 1200, height = 400, res=180)
print (ggseqlogo(l, method='bits'))
dev.off()

#take only the repetitive ones for logo plot
l = list()
ns = c(2386, 2636, 2412, 2977, 16262, 4517)
ptns = c("AC->CA", "CA->AC", "AT->TA", "TA->AT", "AA->TT", "AA->CC")
for (ptn in ptns){
  t = read.table(paste0("~/Downloads/",ptn, "_d1_reponly.tsv"), sep="\t")
  t = as.character(t$V1) #make it a vector
  l = list.append(l,t)
}
names(l) = paste0(ptns, ", n=", ns)
png(paste0("~/Downloads/bits_logo_fig3_final_reponly.png"), width = 1200, height = 400, res=180)
print (ggseqlogo(l, method='bits'))
dev.off()

#d=2, 3, 4, AA->TT / CC->AA
l = list()
ns = c(28882,5990,21533,5204,2629,4110)
names = c()
ptns = c("AA->TT","CC->AA")
ds = c(2,3,4)
for (ptn in ptns){
  for (d in ds){
    t = read.table(paste0("~/Downloads/",ptn, "_d",d,"_reponly.tsv"), sep="\t")
    t = as.character(t$V1) #make it a vector
    l = list.append(l,t)
    names = c(names, paste0(ptn, ", d=",d,", n="))
  }
}
names(l) = paste0(names, ns)
png(paste0("~/Downloads/bits_logo_supp.png"), width = 1200, height = 400, res=180)
print (ggseqlogo(l, method='bits'))
dev.off()





