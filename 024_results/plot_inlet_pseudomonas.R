# vim: set filetype=r
library(ggtree);
library(phangorn);
#library(colorspace);
#library(RColorBrewer);
library(wesanderson);

#tre<-unroot(read.tree("Pseudomonas.fasta.treefile"))
tre<-unroot(read.tree("Pseudomonas_upgma.treefile"))

tre$tip.label <- gsub("Pseudomonas", "P", tre$tip.label) # replace genus by first letter
tre$tip.label <- gsub("\\.\\d{2}$","", tre$tip.label) # remove number with paralog number

toprune = c(
"mendocina", "oryzihabitans", "psychrotolerans", "alcaliphila", "oleovorans", "stutzeri",
"simiae", "orientalis", "frederiksbergensis", "sp000316175", "sp000346225", "sp000349845",
"piscium", "protegens", "ficuserectae", "alkylphenolica", "rhizosphaerae", "cichorii",
"\\.fulva\\_B", "parafulva\\_A", "putida\\_S", "putida\\_J", "putida\\_M",
"citronellolis", "cremoricolorata\\_A", "resinovorans\\_A", "alcaligenes", "soli", "mosselii",
"avellanae", "taiwanensis", "furukawaii", "knackmussii", "putida\\_B", "putida\\_T", "putida\\_E",
"sp000829415", "monteilii\\_A", "marginalis", "synxantha", "lurida", "poae", "monteilii_B\\.0005",
"azotoformans", "fluorescens\\_AA", "fluorescens\\_AN", "fluorescens\\_B",
"entomophila", "viridiflava",
"balearica", "zhaodongensis", "aeruginosa")

remaining = c(
"antarctica",
"brassicacearum",
"corrugata",
"fragi",
"kilonensis",
"mandelii",
"mucidolens",
"sp000931465",
"sp001511755",
"sp001547895",
"sp002966775",
"taetrolens",
"trivialis_B",
"veronii",
"versuta",
"yamanorum"
)

## involved in deep paralogies (i.e. mess the taxon classification depending on copy)
must_remain = c( 
"\\.fulva\\.", "parafulva\\_B",
"brassicacearum",
"cerasi",
"fluorescens",
"fluorescens\\_E",
"fluorescens\\_X",
"fluorescens\\_AR",
"koreensis_B",
"koreensis_C",
"koreensis_D",
"hunanensis",
"syringae",
"syringae_F",
"syringae_M",
"monteilii",
"monteilii_B\\.0003",
"putida",
"putida_H",
"putida_P",
"putida_Q",
"chlororaphis",
"sp001655595",
"sp001655615"
)

for (i in 1:length(toprune)) { 
  tre <- drop.tip (tre, grep (toprune[i], tre$tip.label));
}

#x=quantile(tre$edge.length, p=0.99)[[1]];
#for (j in 1:length(tre$edge.length)) {
#  if (tre$edge.length[j] > 2 * x) tre$edge.length[j] = 2 * x
#}
tre <- midpoint (tre) ## after outlier truncation 

# nested gsub removes prefix and sufix: A_B_C --> B
groupInfo <- split(tre$tip.label, gsub("\\.\\w+","", gsub("^\\w*?\\.","", tre$tip.label)) )
grouptre <- groupOTU(tre,groupInfo)
color_palette <- wes_palette("FantasticFox1", length(groupInfo) + 5, type="continuous"); # FantasticFox1 Cavalcanti1

p <- ggtree(grouptre, layout="fan", open.angle=171, size=0.8, aes(color=group)) # angle=271 
p <- p + geom_tiplab2(align=TRUE, linetype=1, linesize=0.15, fontface = "bold", size=1.2)  
p <- p + ggplot2::xlim(-0.05,NA)
p <- p + ggplot2::scale_color_manual(values=color_palette) 
ggplot2::ggsave ("Pseudomonas_inlet.png", width = 16, height= 16, bg = "white");
quit();

#ggtree(mytre, aes(color=group), layout='circular') + geom_tiplab(size=1, aes(angle=angle))
## for arbitrary colorspaces:
# geom_balance and geom_hilight  to create rectangles around nodes
#color_palette  <- rainbow_hcl(length(groupInfo) + 1) 
#color_palette  <- colorRampPalette(brewer.pal(7, "Paired"))(length(groupInfo) + 1);
# p <- ggtree(grouptre, layout="circular", aes(color=group), size=0.9) 
# p <- p + geom_tiplab(linetype=1, offset=0, size=1.2) 
# p <- ggtree(grouptre, aes(color=group), branch.length="none") 
# p <- p + geom_tiplab(align=TRUE, aes(angle=angle), linetype=1, linesize=0.1, size=1) # tiplab2() adjusts angle better 
# p <- ggtree(grouptre, layout="fan", open.angle=5, aes(color=group)) 
