# vim: set filetype=r
library(ape);
library(RColorBrewer);

png ("concat_tree%02d.png", width = 6000, height= 10000, bg = "white");
par (tcl=-0.15, mgp=c(2.1,0.2,0), mar=c(1,1,0.1,0.1));
par (font=2, cex=2, col.main="gray10", cex.axis=1.5, cex.main=5, cex.lab=1.2);
layout(matrix(1:2, 1, 2, byrow = TRUE))


minclad <- c("Campylobacter", "Clostridium", "Enterococcus", "Escherichia", "Helicobacter", "Leptospira", "Listeria", 
             "Mycobacterium", "Neisseria", "Pseudomonas", "Salmonella", "Staphylococcus", "Streptococcus");
# minclad <- c("Enterococcus","Escherichia", "Pseudomonas","Salmonella","Staphylococcus","Streptococcus");
mic <- colorRampPalette(brewer.pal(13, "Dark2"))(length(minclad));

treefile <- list.files(pattern="*.treefile")
for (i in 1:length(treefile)) {
    tre<-read.tree(treefile[i])
    len<-length(tre$tip.label);
    col1 <- rep("black",len);
    for (k in 1:len) { # notice that we treat Enterococcus and Enterococcus_B as different genera/colours...
        for (j in 1:length(minclad)) { # look for anything at the beginning except ".", followed by "minclad[].".
            if ( grepl(paste("^",minclad[j],".",sep=""),tre$tip.label[k]) ) col1[k] <- mic[j]; 
        }
        x = unlist(strsplit(tre$tip.label[k], "."));
        #tre$tip.label[k] = paste(x[1],x[2]); # genus + sp 
    }
    my_cex = ifelse (len < 300, 3, 0.9);
    plot (tre, cex=my_cex, node.depth=1, font=2, edge.width=12, tip.color=col1, no.margin=T, label.offset=0.001, direction="rightwards");
    x = unlist(strsplit(treefile[i],'[.]treefile'))  # split is regexp, so dot is special
    title(main=x[1], line = -3);
    legend("bottomleft", legend=minclad, bty="n", fill=mic, border=mic, cex=5, inset=0.1);
}



dev.off(); quit();

groupInfo <- split(tre$tip.label, gsub(".\\w+","", tre$tip.label))
mytre <- groupOTU(tre,groupInfo)
ggtree(mytre, aes(color=group), layout='circular') + geom_tiplab(size=1, aes(angle=angle))
## for arbitrary colorspaces:
library("colorspace")
ggtree(tree, aes(color=group, linetype=group)) + geom_tiplab() + scale_color_manual(values=c("black", rainbow_hcl(4))) + theme(legend.position="right")
# geom_balance and geom_hilight  to create rectangles around nodes
