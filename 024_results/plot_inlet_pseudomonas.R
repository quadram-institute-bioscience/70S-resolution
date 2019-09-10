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
"ficuserectae", "fluorescens_AA", "fluorescens_E", "fluorescens_B", "frederiksbergensis_B", "frederiksbergensis_E", 
"koreensis_B", "koreensis_C", "koreensis_D", "lurida", "mandelii", "marginalis", "mosselii",
"mucidolens", "orientalis_A","proteolytica", "piscium", "putida_S", "simiae", "parafulva_A", "poae",
"protegens", "trivialis_B", "veronii", "silesiensis", "knackmussii", "psychrotolerans",
"cichorii", "rhizosphaerae", "alkylphenolica", "citronellolis", "cremoricolorata_A", "kilonensis",
"sp000346225", "sp000931465", "sp001511755", "sp000829415", "sp001547895", "sp001655615",
"sp001708505", "sp002753995", "sp002966775", "sp002874965", "sp003151075", "sp000349845",
"sp002163625",
"aeruginosa_A", "oryzihabitans_E", "alcaligenes","balearica", "furukawaii",
"alcaliphila_B", "oleovorans", "fulva_B", "mendocina", "mendocina_A", "mendocina_C",
"toyotomiensis", "taiwanensis", "taetrolens", "viridiflava", "yamanorum", "brassicacearum_A"
)

remaining = c(
"fragi",
"fulva",
"hunanensis",
"monteilii",
"monteilii_A",
"monteilii_B",
"parafulva_B",
"putida",
"putida_B",
"putida_E",
"putida_H",
"putida_J",
"putida_M",
"putida_P",
"putida_Q",
"putida_T",
"resinovorans_A",
"sp000316175",
"sp002025205",
"sp002056295",
"sp003205815",
"stutzeri",
"stutzeri_A",
"stutzeri_AE",
"stutzeri_AG",
"stutzeri_AJ",
"stutzeri_D",
"stutzeri_U",
"versuta",
"zhaodongensis"
)

## involved in deep paralogies (i.e. mess the taxon classification depending on copy)
must_remain = c(  
"cerasi",
"brassicacearum_C",
"synxantha",
"syringae",
"syringae_F",
"syringae_M",
"lactis",
"fluorescens",
"fluorescens_AN",
"fluorescens_AR",
"fluorescens_X",
"sp001655595",
"chlororaphis",
"chlororaphis_E"
)

for (i in 1:length(toprune)) { 
  tre <- drop.tip (tre, grep (toprune[i], tre$tip.label));
}

x=quantile(tre$edge.length, p=0.99)[[1]];
for (j in 1:length(tre$edge.length)) {
  if (tre$edge.length[j] > 2 * x) tre$edge.length[j] = 2 * x
}
tre <- midpoint (tre) ## after outlier truncation 

# nested gsub removes prefix and sufix: A_B_C --> B
groupInfo <- split(tre$tip.label, gsub("\\.\\w+","", gsub("^\\w*?\\.","", tre$tip.label)) )
grouptre <- groupOTU(tre,groupInfo)
#color_palette <- wes_palette("FantasticFox1", length(groupInfo) + 5, type="continuous"); # FantasticFox1 Cavalcanti1
color_palette  <- colorRampPalette(brewer.pal(12, "Paired"))(length(groupInfo) + 1);

p <- ggtree(grouptre, layout="fan", open.angle=100, size=0.8, aes(color=group)) # angle=271 
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
