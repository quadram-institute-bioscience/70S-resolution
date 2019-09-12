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
"avellanae", "antarctica\\_A", "azotoformans\\_A", "azotoformans\\_B", "corrugata", "entomophila",
"ficuserectae", "fluorescens\\_AA", "fluorescens\\_E", "fluorescens\\_B", "frederiksbergensis\\_B", "frederiksbergensis\\_E", 
"koreensis_B", "koreensis_C", "koreensis_D", "lurida", "mandelii", "marginalis", "mosselii",
"mucidolens", "orientalis\\_A","proteolytica", "piscium", "putida\\_S", "simiae", "parafulva\\_A", "\\.poae\\.",
"protegens", "trivialis\\_B", "veronii", "silesiensis", "knackmussii", "psychrotolerans",
"cichorii", "rhizosphaerae", "alkylphenolica", "citronellolis", "cremoricolorata\\_A", "kilonensis",
"sp000346225", "sp000931465", "sp001511755", "sp000829415", "sp001547895", "sp001655615",
"sp001708505", "sp002753995", "sp002966775", "sp002874965", "sp003151075", "sp000349845", "sp002163625",
"aeruginosa\\_A", "oryzihabitans\\_E", "alcaligenes","balearica", "furukawaii",
"alcaliphila\\_B", "oleovorans", "\\.fulva\\_B", "mendocina", "mendocina\\_A", "mendocina\\_C",
"fluorescens\\.","fluorescens\\_AN","fluorescens\\_AR","fluorescens\\_X", "chlororaphis","chlororaphis\\_E", "brassicacearum\\_C", 
"sp000316175","cerasi","lactis","versuta","fragi","sp001655595","syringae\\.","syringae\\_F","syringae\\_M",
"monteilii\\_B\\.0005","monteilii\\.","monteilii\\_A","synxantha","fulva\\.001","putida\\_P\\.001",
"putida\\_J","putida\\_M","putida\\_B","putida\\_E","putida\\_T","sp002056295", "sp003205815",
"stutzeri\\_A","stutzeri\\_AG","stutzeri\\_D","stutzeri\\_U","stutzeri\\_AE","stutzeri\\_AJ","resinovorans\\_A","sp002025205",
"zhaodongensis","stutzeri","toyotomiensis", "taiwanensis", "taetrolens", "viridiflava", "yamanorum", "brassicacearum\\_A"
)

remaining = c(
"hunanensis",
"monteilii",
"monteilii\\_B",
"putida",
"putida\\_H",
"putida\\_P",
"putida\\_Q"
)

#"fulva","parafulva\\_B",

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
color_palette <- wes_palette("Cavalcanti1", length(groupInfo) + 1, type="continuous"); # FantasticFox1 Cavalcanti1

p <- ggtree(grouptre, layout="fan", open.angle=292, size=0.8, aes(color=group)) # angle=271 
p <- p + geom_tiplab2(align=TRUE, linetype=1, linesize=0.15, fontface = "bold", size=3.2)
p <- p + ggplot2::xlim(-0.05,0.1)
p <- p + ggplot2::scale_color_manual(values=color_palette) 
ggplot2::ggsave ("Pseudomonas_inlet_2.png", width = 16, height= 16, bg = "white");
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
