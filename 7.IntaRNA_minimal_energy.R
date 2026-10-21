# Minimal energy profiles
#pdf(file="Minimal_energy_ADH4.pdf")
d <- read.table("MYPROFILEFILE_SPARC.csv", header=T, sep=";");
plot( d[,1], d[,3], xlab="SPARC", ylab="minimal energy", type="l", col="blue", lwd=2)
abline(h=0, col="red", lty=2, lwd=2)
#dev.off()

#pdf(file="Minimal_energy_lnc-ADH4-3-2.pdf")
d <- read.table("MYPROFILEFILE_lnc_SPARC.csv", header=T, sep=";");
plot( d[,1], d[,3], xlab="lnc-SPARC-2:3", ylab="minimal energy", type="l", col="blue", lwd=2)
abline(h=0, col="red", lty=2, lwd=2)
#dev.off()

# Minimal energy for all inter molecular index pairs
library(tidyverse)
#pdf(file=" intermolecular_index_pairs_lnc-ADH4-3-2.pdf")
read_delim("MYPAIRMINE_SPARC.csv", delim=";") %>%
  replace(is.na(.), 0) %>%
  rename( target = minE ) %>%
  pivot_longer( cols=-target, names_to = "query", values_to = "E" ) %>%
  separate( query, into = c("qSeq","qIdx"), sep = "_") %>%
  separate( target, into = c("tSeq","tIdx"), sep = "_") %>%
  mutate(across(c(tIdx,qIdx,E),as.numeric)) %>%
  mutate(E = ifelse( E>0, 0, E)) %>%
  ggplot(aes(x=tIdx,y=qIdx, fill=E)) +
  geom_tile() + 
  xlab("SPARC") + ylab("lnc-SPARC-2:3")+
  scale_fill_gradient(low="red",high="yellow")
