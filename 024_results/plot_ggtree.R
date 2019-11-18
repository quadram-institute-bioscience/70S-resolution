# vim: set filetype=r
library(ggtree);
library(phangorn);
library(colorspace);
library(RColorBrewer);

## leaf names are Pseudomonas_E.hunanensis.002953115.04 ( genus . sp . accession . paralognumber )
treefile <- list.files(pattern="fasta.treefile")
#treefile <- list.files(pattern="_upgma.treefile")
for (i in 1:length(treefile)) {
    tre<-unroot(read.tree(treefile[i]))
    x=quantile(tre$edge.length, p=0.999)[[1]];
    for (j in 1:length(tre$edge.length)) {
        if (tre$edge.length[j] > 1 * x) tre$edge.length[j] = 1 * x
    }
    tre <- midpoint (tre) ## after outlier truncation  (not needed when plotting upgma trees) 

    tre$tip.label <- gsub("Pseudomonas", "P", tre$tip.label)
    tre$tip.label <- gsub("Staphylococcus", "S", tre$tip.label) # replace genus by first letter
    tre$tip.label <- gsub("\\.\\d{2}$","", tre$tip.label) # remove number with paralog number

    # remove (tre, which (tre$tip.label %in% "900240075"))
    # nested gsub removes prefix and sufix: A_B_C --> B
    groupInfo <- split(tre$tip.label, gsub("\\.\\w+","", gsub("^\\w*?\\.","", tre$tip.label)) )
    grouptre <- groupOTU(tre,groupInfo)
    color_palette  <- colorRampPalette(brewer.pal(12, "Paired"))(length(groupInfo) + 1);
    #color_palette  <- rainbow_hcl(length(groupInfo) + 1) 
    fname = unlist(strsplit(treefile[i],'[.]fasta[.]treefile'))[1]  # split is regexp, so dot is special
    # fname = unlist(strsplit(treefile[i],'_upgma[.]treefile'))[1]  # split is regexp, so dot is special

    tscaling = 1.; blen = 0.55; 
    if (length(tre$tip.label) < 200) { tscaling = 2.4; blen = 0.22; }

    #p <- ggtree(grouptre, layout="circular", aes(color=group), size=0.9) 
    p <- ggtree(grouptre, layout="fan", open.angle=70, size=1. * tscaling, aes(color=group)) 
    p <- p + geom_tiplab2(align=TRUE,linetype=1, linesize=0.15, size=1.4 * tscaling)  
    p <- p + ggplot2::scale_color_manual(values=color_palette) 
    p + ggplot2::ggtitle(fname)
    ggplot2::ggsave (paste(fname,"_ggcircular.png",sep=""), width = 20, height= 20, bg = "white");


    tre$edge.length <- log(tre$edge.length + 1.)
    color_palette  <- colorRampPalette(brewer.pal(7, "Dark2"))(length(groupInfo) + 1); ## or Dark2 or Set1
    p <- ggtree(grouptre, aes(color=group), size = 1. * tscaling) 
    p <- p + geom_tiplab(align=TRUE, linetype=1, linesize=0.2, size=0.9 * tscaling) 
    p <- p + ggplot2::scale_color_manual(values=color_palette) 
    p <- p + theme_tree2() # adds bottom x scale 
    p <- p + ggplot2::xlim(0, blen);
    p + ggplot2::ggtitle(fname)

    ggplot2::ggsave (paste(fname,"_ggtree.png",sep=""), width = 10, height= 20, bg = "white");
}
quit();

#ggtree(mytre, aes(color=group), layout='circular') + geom_tiplab(size=1, aes(angle=angle))
## for arbitrary colorspaces:
# geom_balance and geom_hilight  to create rectangles around nodes

# p <- p + geom_tiplab(linetype=1, offset=0, size=1.2) 
# p <- ggtree(grouptre, aes(color=group), branch.length="none") 
# p <- p + geom_tiplab(align=TRUE, aes(angle=angle), linetype=1, linesize=0.1, size=1) # tiplab2() adjusts angle better 
# p <- ggtree(grouptre, layout="fan", open.angle=5, aes(color=group)) 
# p <- p + theme(legend.position="right")
