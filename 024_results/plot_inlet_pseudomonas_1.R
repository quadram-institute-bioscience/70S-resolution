# vim: set filetype=r
library(ggtree);
library(phangorn);
library(colorspace);
library(RColorBrewer);
library(wesanderson);

tre<-unroot(read.tree("Pseudomonas.fasta.treefile"))

tre$tip.label <- gsub("Pseudomonas", "P", tre$tip.label) # replace genus by first letter
tre$tip.label <- gsub("\\.\\d{2}$","", tre$tip.label) # remove number with paralog number

toprune = c(
"avellanae", "antarctica_A", "azotoformans_A", "azotoformans_B", "corrugata", "entomophila",
"ficuserectae", "fluorescens\\_AA", "fluorescens\\_E", "fluorescens\\_B", "frederiksbergensis\\_B", "frederiksbergensis\\_E", 
"koreensis\\_B", "koreensis\\_C", "koreensis\\_D", "lurida", "mandelii", "marginalis", "mosselii",
"mucidolens", "orientalis_A","proteolytica", "piscium", "putida_S", "simiae", "parafulva\\_A", "poae",
"protegens", "trivialis_B", "veronii", "silesiensis", "knackmussii", "psychrotolerans",
"cichorii", "rhizosphaerae", "alkylphenolica", "citronellolis", "cremoricolorata\\_A", "kilonensis",
"sp000346225", "sp000931465", "sp001511755", "sp000829415", "sp001547895", "sp001655615",
"sp001708505", "sp002753995", "sp002966775", "sp002874965", "sp003151075", "sp000349845",
"sp002163625", "aeruginosa_A", "oryzihabitans_E", "alcaligenes","balearica", "furukawaii",
"alcaliphila\\_B", "oleovorans", "\\.fulva\\_B", "mendocina", "mendocina\\_A", "mendocina\\_C",
"toyotomiensis", "taiwanensis", "taetrolens", "viridiflava", "yamanorum", "brassicacearum\\_A",

"syringae\\.003", "syringae\\.0009","hunanensis","fragi","\\.fulva\\.","\\.parafulva\\_B","brassicacearum\\_C",
"\\.fluorescens\\.","sp002025205","sp002056295","sp003205815","sp000316175","cerasi\\.001","chlororaphis\\.",
"monteilii","monteilii\\_A","monteilii\\_B","versuta","zhaodongensis","resinovorans_A","syringae\\_F",
"putida","putida\\_B", "putida\\_E","putida\\_H","putida\\_J", "putida\\_M", "putida\\_P", "putida\\_Q", "putida\\_T",
"stutzeri", "stutzeri\\_A", "stutzeri\\_AE", "stutzeri\\_AG", "stutzeri\\_AJ", "stutzeri\\_D", "stutzeri\\_U"
)


## <<var not used>> involved in deep paralogies (i.e. mess the taxon classification depending on copy)
must_remain = c(  
"synxantha",
"syringae",
"syringae_M",
"lactis",
"fluorescens_AN",
"fluorescens_AR",
"fluorescens_X",
"sp001655595",
"chlororaphis_E"
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
color_palette <- wes_palette("Cavalcanti1", length(groupInfo) + 2, type="continuous"); # FantasticFox1 Cavalcanti1

p <- ggtree(grouptre, layout="fan", open.angle=280, size=1., aes(color=group)) # angle=271 
p <- p + geom_tiplab2(align=TRUE, linetype=1, linesize=0.15, fontface = "bold", size=3.)  
p <- p + ggplot2::xlim(-0.09,0.26)
p <- p + ggplot2::scale_color_manual(values=color_palette) 
ggplot2::ggsave ("Pseudomonas_inlet_1.png", width = 16, height= 16, bg = "white");
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
