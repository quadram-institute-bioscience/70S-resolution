# vim: set filetype=r
library(ggtree);
library(phangorn);
library(colorspace);
library(RColorBrewer);

treefile <- list.files(pattern="fasta.treefile")
for (i in 1:length(treefile)) {
    tre<-unroot(read.tree(treefile[i]))
    x=quantile(tre$edge.length, p=0.99)[[1]];
    for (j in 1:length(tre$edge.length)) {
        if (tre$edge.length[j] > 1 * x) tre$edge.length[j] = 1 * x
    }
    tre <- midpoint (tre) ## after outlier truncation 
    tre$tip.label <- gsub("\\.\\d+$","", tre$tip.label) # remove accession number
    #groupInfo <- split(tre$tip.label, gsub("_\\w+","", gsub("^\\w*?_","", tre$tip.label)) )# nested gsub removes prefix and sufix: A_B_C --> B
    groupInfo <- split(tre$tip.label, gsub("\\.\\w+","", tre$tip.label) ) # color by genus
    grouptre <- groupOTU(tre,groupInfo)
    color_palette  <- colorRampPalette(brewer.pal(12, "Paired"))(length(groupInfo) + 1);
    fname = unlist(strsplit(treefile[i],'[.]fasta[.]treefile'))[1]  # split is regexp, so dot is special

    p <- ggtree(grouptre, layout="circular", aes(color=group), size=1.) 
    p <- p + geom_tiplab2(align=TRUE,linetype=1, linesize=0.15, size=2.)  
    p <- p + ggplot2::scale_color_manual(values=color_palette) 
    p <- p + ggplot2::xlim(-0.01,NA)
    p + ggplot2::ggtitle(fname)
    ggplot2::ggsave (paste(fname,"_ggcircular.png",sep=""), width = 20, height= 20, bg = "transparent");

    tre$edge.length <- log(tre$edge.length + 1.)
    color_palette  <- colorRampPalette(brewer.pal(7, "Dark2"))(length(groupInfo) + 1); ## or Dark2 or Set1
    p <- ggtree(grouptre, aes(color=group)) 
    p <- p + geom_tiplab(align=FALSE, linetype=1, linesize=0.2, size=1.1) 
    p <- p + ggplot2::scale_color_manual(values=color_palette) 
    p <- p + theme(legend.position="right")
    p + ggplot2::ggtitle(fname)
    ggplot2::ggsave (paste(fname,"_ggtree.png",sep=""), width = 8, height= 20, bg = "white");
}
quit();

# geom_balance and geom_hilight  to create rectangles around nodes
# color_palette  <- rainbow_hcl(length(groupInfo) + 1) 
# p <- ggtree(grouptre, layout="fan", open.angle=70, size=0.8, aes(color=group)) 
# p <- ggtree(grouptre, aes(color=group), branch.length="none") 
# p <- p + geom_tiplab(align=TRUE, aes(angle=angle), linetype=1, linesize=0.1, size=1) # tiplab2() adjusts angle better 
# p <- p + geom_tiplab(linetype=1, offset=0, size=1.2) 
# p <- p + theme_tree2() # adds bottom x scale 
