
plot_for_manu_tropism2 <- function(sample1, ObjectPX5cds, group.name, date1 ="0812_2021") {
obj   <- subset(ObjectPX5cds, subset = Finalname %in% sample1)
obj@meta.data$Finalname <- factor(obj@meta.data$Finalname , levels = sample1)

data_df <- data.frame(obj@reductions$umap@cell.embeddings)
data_df$pseudotime <- obj@meta.data$pseudotime
data_df$Finalname <- obj@meta.data$Finalname
data_df$Cluster <- obj@meta.data$seurat_clusters
data_df$Phase <- obj@meta.data$Phase
head(data_df)

p1 <- ggplot(data_df, aes(UMAP_1, UMAP_2)) + 
       geom_point(aes(colour = Cluster), size = 1) + #2
	   #scale_color_jcolors_contin("pal3", bias = 1.75) + 
	   #facet_wrap( ~ Finalname, nrow=1) +  
	   theme_classic() +
	   theme(axis.text=element_text(size=10), axis.title=element_text(size=14))   +
	   #theme(strip.text.x = element_text(size = 20)) +
	   #theme(strip.background = element_blank()) +
	   theme(legend.title = element_text(size = 20)) +
	   theme(legend.text = element_text(size=15)) + guides(colour = guide_legend(override.aes = list(size=5))) 

png(paste0("figure/",date1, ".monocle3.", group.name, "-clusters1.png"), width = 520, height =460)
print(p1)
dev.off()

p2 <- ggplot(data_df, aes(UMAP_1, UMAP_2)) + 
       geom_point(aes(colour = Phase), size = 1) +   #2
	   #scale_color_jcolors_contin("pal3", bias = 1.75) + 
	   #facet_wrap( ~ Finalname, nrow=1) +  
	   theme_classic() +
	   theme(axis.text=element_text(size=10), axis.title=element_text(size=14))   +
	   #theme(strip.text.x = element_text(size = 20)) +
	   #theme(strip.background = element_blank()) +
	   theme(legend.title = element_text(size = 20)) +
	   theme(legend.text = element_text(size=15)) + guides(colour = guide_legend(override.aes = list(size=5)))

png(paste0("figure/",date1, ".monocle3.", group.name, "-Phase1.png"), width = 520, height =460)
print(p2)
dev.off()

p3 <- ggplot(data_df, aes(UMAP_1, UMAP_2)) + 
       geom_point(aes(colour =Cluster), size = 1) + 
	   #scale_color_jcolors_contin("pal3", bias = 1.75) + 
	   facet_wrap( ~ Finalname, nrow=1) +  
	   theme_classic() +
	   theme(axis.text=element_text(size=10), axis.title=element_text(size=14))   +
	   #theme(strip.text.x =  element_blank()) +   
	   theme(strip.text.x = element_text(size = 15)) +
	   theme(strip.background = element_blank()) +
	   theme(legend.title = element_text(size = 18)) +
	   theme(legend.text = element_text(size=15)) + 
	   guides(colour = guide_legend(override.aes = list(size=5)))
png(paste0("figure/",date1, ".monocle3.", group.name, "-clusters.png"), width = 800, height =200)
print(p3)
dev.off()

p4 <- ggplot(data_df, aes(UMAP_1, UMAP_2)) + 
       geom_point(aes(colour =Phase), size = 1) + 
	   #scale_color_jcolors_contin("pal3", bias = 1.75) + 
	   facet_wrap( ~ Finalname, nrow=1) +  
	   theme_classic() +
	   theme(axis.text=element_text(size=10), axis.title=element_text(size=14))   +
	   #theme(strip.text.x =  element_blank()) +  
	   theme(strip.text.x = element_text(size = 15)) +
	   theme(strip.background = element_blank()) +
	   theme(legend.title = element_text(size = 18)) +
	   theme(legend.text = element_text(size=15)) + 
	   guides(colour = guide_legend(override.aes = list(size=5)))

png(paste0("figure/",date1, ".monocle3.", group.name, "-Phase.png"), width = 800, height =200)
print(p4)
dev.off()


######### Phase_barplot(obj=obj, sample1= sample1, group.name= group.name, date1 ="0812_2021")

ident_inf <- unique(obj@meta.data$seurat_clusters)  #unique(Idents(obj))  # unique(obj@meta.data$clusters)
Source_inf <- unique(obj@meta.data$Source)
Finalname_inf <- unique(obj@meta.data$Finalname)
Pt_inf <- unique(obj@meta.data$Pt)
model_inf <- unique(obj@meta.data$model)
Phase_inf <- unique(obj@meta.data$Phase)

alist <- c()
for (i in 1:length(Finalname_inf)) {		 
alist <- c(alist, subset(obj, subset = Finalname %in% Finalname_inf[i]))
}

dataset <- alist[[1]]
clist <- c()

for (j in ident_inf) {		 
clist <- c(clist, nrow(dataset@meta.data[dataset@meta.data$seurat_clusters %in% j,])/ncol(dataset))
}

 data1 <- t(data.frame(clist))
 colnames(data1) <- ident_inf
 rownames(data1) <- Finalname_inf[1]
 data2 <- data1

for (i in 2:length(Finalname_inf)) {	
	 #print(i)
	 dataset <- alist[[i]]
	 #print(dataset)
clist <- c()	
for (j in ident_inf) {		 
clist <- c(clist, nrow(dataset@meta.data[dataset@meta.data$seurat_clusters %in% j,])/ncol(dataset))
}

 data1 <- t(data.frame(clist))
 colnames(data1) <- ident_inf
 rownames(data1) <- Finalname_inf[i]
	data2 <- rbind(data2, data1)	
	 }
	 
data2 <-  data2[sample1,sort(colnames(data2))]
print(round(data2 ,4) )

# write.table(round(data2 ,2), file= "0319.fastMNN.txt")
	 library("reshape2")
	 data3.Finalname <- melt(data2)
	 data3.Finalname$Var2 <- as.character(data3.Finalname$Var2)

colnames(data3.Finalname) <- c("Name", "Cluster", "Percentage")
data3.Finalname$Percentage <- data3.Finalname$Percentage  * 100
#print(data3.Finalname)

#data3.Finalname$Name <- factor(data3.Finalname$Name, levels = sample1)

p5 <- ggplot(data3.Finalname, aes(Name, Percentage, fill = Cluster)) +
      geom_bar(position="stack", stat="identity") +
      scale_x_discrete(guide = guide_axis(angle = 45))  + 
      theme_classic()+
	  theme(axis.text=element_text(size=10), axis.title=element_text(size=14))   +
      theme(legend.title = element_text(size =12 )) +
	  theme(legend.text = element_text(size=8)) +
	  theme(axis.title.x=element_blank()) + 
	  guides(colour = guide_legend(override.aes = list(size=5))) + theme(legend.key.size = unit(5, 'mm'))

png(paste0("./figure/", date1, "_", group.name,  "_Finalname.clusters.barplot.png"), width = 520, height = 460)
print( p5 )
dev.off()




############cluster_barplot(obj=obj, sample1= sample1, group.name= group.name, date1 ="0812_2021")


alist <- c()
for (i in 1:length(Finalname_inf)) {		 
alist <- c(alist, subset(obj, subset = Finalname %in% Finalname_inf[i]))
}

dataset <- alist[[1]]
clist <- c()

for (j in Phase_inf) {		 
clist <- c(clist, nrow(dataset@meta.data[dataset@meta.data$Phase %in% j,])/ncol(dataset))
}

 data1 <- t(data.frame(clist))
 colnames(data1) <- Phase_inf
 rownames(data1) <- Finalname_inf[1]
 data2 <- data1

for (i in 2:length(Finalname_inf)) {	
	 #print(i)
	 dataset <- alist[[i]]
	 #print(dataset)
clist <- c()	
for (j in Phase_inf) {		 
clist <- c(clist, nrow(dataset@meta.data[dataset@meta.data$Phase %in% j,])/ncol(dataset))
}

 data1 <- t(data.frame(clist))
 colnames(data1) <- Phase_inf
 rownames(data1) <- Finalname_inf[i]
	data2 <- rbind(data2, data1)	
	 }
	 
data2 <-  data2[sample1,sort(colnames(data2))]
print(round(data2 ,4) )

# write.table(round(data2 ,2), file= "0319.fastMNN.txt")
	 library("reshape2")
	 data3.Finalname <- melt(data2)
	 data3.Finalname$Var2 <- as.character(data3.Finalname$Var2)

colnames(data3.Finalname) <- c("Name", "Phase", "Percentage")
data3.Finalname$Percentage <- data3.Finalname$Percentage  * 100
#print(data3.Finalname)

#data3.Finalname$Name <- factor(data3.Finalname$Name, levels = sample1)

p6 <- ggplot(data3.Finalname, aes(Name, Percentage, fill = Phase))+
      geom_bar(position="stack", stat="identity") +
      scale_x_discrete(guide = guide_axis(angle = 45))  + 
      theme_classic()+
	  theme(axis.text=element_text(size=10), axis.title=element_text(size=14))   +
      theme(legend.title = element_text(size = 10)) +
	  theme(legend.text = element_text(size=8)) +
	  theme(axis.title.x=element_blank()) + 
	  guides(colour = guide_legend(override.aes = list(size=5)))
	  
png(paste0("./figure/", date1, "_", group.name,  "_Finalname.Phase.barplot.png"), width = 520, height = 460)
print(p6)
dev.off()

p35 <- ggarrange(p1 + theme(legend.position = "none") + ggtitle("Cluster") + theme(plot.title = element_text(color="red", size=18 )) + theme( plot.title = element_text(hjust = 0.5)) ,
                 p3 + theme(legend.position = "none"), 
                 p5,  widths  = c(1, 4, 1.2),
          ncol = 3, nrow = 1)

p46 <- ggarrange(p2 + theme(legend.position = "none") + ggtitle("Phase") + theme(plot.title = element_text(color="red", size=18 )) + theme( plot.title = element_text(hjust = 0.5)),
                 p4 + theme(legend.position = "none"),
                 p6,  widths  = c(1, 4, 1.2),
          ncol = 3, nrow = 1)
	  
	  
png(paste0("./figure/", date1, "_", group.name,  "_Finalname.1-6.png"), width =1200, height =400)
print(
ggarrange(p35, p46, heights = c(1, 1),
          ncol = 1, nrow = 2)
)
dev.off()
}

