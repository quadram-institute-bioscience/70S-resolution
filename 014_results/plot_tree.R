# vim: set filetype=r
library(ape);
library(RColorBrewer);

png ("operon%02d.png", width = 6000, height= 12000, bg = "white");
par (tcl=-0.15, mgp=c(2.1,0.2,0), mar=c(1,1,0.1,0.1));
par (font=2, cex=2, col.main="gray10", cex.axis=1.5, cex.main=5, cex.lab=1.2);
#layout(matrix(1:2, 1, 2, byrow = TRUE))

treefile <- list.files(pattern="fasta.treefile")
for (i in 1:length(treefile)) {
    tre<-read.tree(treefile[i])
    len<-length(tre$tip.label);
    name_table = matrix(ncol=2, nrow=len)
    for (k in 1:len) {
				x = unlist(strsplit(tre$tip.label[k], "_"));
        name_table[k,1] = x[2]
        name_table[k,2] = paste(x[2], x[3], sep = "_")
    }
    # each tip label is: genus_species_strain_paralog
#   sp_species <- name_table[,2]; names (sp_species) <- tre$tip.label; # gives species for each leaf label
#   sp_strain  <- name_table[,3]; names (sp_strain)  <- tre$tip.label; # gives strain (sample id) for each leaf
    majclad <- unique(name_table[,1]); # species
    minclad <- unique(name_table[,2]); # sample = species_strain
    mac <- colorRampPalette(brewer.pal(12, "Paired"))(length(majclad));
    mic <- colorRampPalette(brewer.pal(8, "Dark2"))(length(minclad));
    col_mac <- rep("black",len);
    col_mic <- rep("black",len);

    for (k in 1:len) {
        for (j in 1:length(majclad)) { 
            if ( grepl(paste(majclad[j],"_",sep=""),tre$tip.label[k]) ) col_mac[k] <- mac[j];
        }
        for (j in 1:length(minclad)) { 
            if ( grepl(paste(minclad[j],"_",sep=""),tre$tip.label[k]) ) col_mic[k] <- mic[j];
        }
    }
    my_cex = ifelse (len < 500, 1.2, .9);
    plot (tre, cex=my_cex, node.depth=1, font=2, edge.width=12, tip.color=col_mac, no.margin=T, label.offset=0.001, direction="rightwards");
    x = unlist(strsplit(treefile[i],'[.]treefile'))  # split is regexp, so dot is special
    title(main=x[1], line = -3);
    if ( i == 1) { legend("bottomleft", legend=majclad, bty="n", fill=mac, border=mac, cex=5, inset=0.1);
    } else { legend("topright", legend=majclad, bty="n", fill=mac, border=mac, cex=5, inset=0.1); }
}
dev.off(); quit();

