###Written by: Yanlei Feng & Susann Wicke
###CITE: Feng & Wicke (2022) First mitochondrial genomes of true ferns allow modeling the mitogenomic inflation syndrome across all land plant lineages. bioarxiv.

a<-read.table("mitogenome_features.tsv", sep='\t', header=T)

library( ggplot2 )
library( gridExtra )
library( ggpubr )

#violin
pdf("violin.pdf")
point_size=0.45
point_color='turquoise2'
grid_line=0.5
size<-ggplot(a, aes(x=Clade2, y=round(size/1000)))+geom_violin()+coord_flip()+theme(axis.ticks.x = element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position = "none",panel.grid.minor.x = element_line(size =grid_line),panel.grid.minor.y = element_line(size =grid_line),panel.grid.major.x = element_line(size =grid_line),panel.grid.major.y = element_line(size =grid_line))+scale_y_continuous(limits=c(0,2000),breaks=c(0,600,1200,1800))+geom_point(color=point_color,size=point_size)+stat_summary(fun.y=median,geom="point",size=1,aes(group=1,color='red'))
gc<-ggplot(a, aes(x=Clade2, y=GC))+geom_violin()+coord_flip()+theme(axis.ticks.x = element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position = "none",panel.grid.minor.x = element_line(size =grid_line),panel.grid.minor.y = element_line(size =grid_line),panel.grid.major.x = element_line(size =grid_line),panel.grid.major.y = element_line(size =grid_line))+scale_y_continuous(limits=c(25,55),breaks=c(30,40,50))+geom_point(color=point_color,size=point_size)+stat_summary(fun.y=median,geom="point",size=1,aes(group=1,color='red'))
genes<-ggplot(a, aes(x=Clade2, y=Genes))+geom_violin()+coord_flip()+theme(axis.ticks.x = element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position = "none",panel.grid.minor.x = element_line(size =grid_line),panel.grid.minor.y = element_line(size =grid_line),panel.grid.major.x = element_line(size =grid_line),panel.grid.major.y = element_line(size =grid_line))+scale_y_continuous(limits=c(17,43),breaks=c(20,30,40))+geom_point(color=point_color,size=point_size)+stat_summary(fun.y=median,geom="point",size=1,aes(group=1,color='red'))
# ggarrange( size, gc, genes, nrow=1, ncol=3 )
intron<-ggplot(a, aes(x=Clade2, y=Intron))+geom_violin()+coord_flip()+theme(axis.ticks.x = element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position = "none",panel.grid.minor.x = element_line(size =grid_line),panel.grid.minor.y = element_line(size =grid_line),panel.grid.major.x = element_line(size =grid_line),panel.grid.major.y = element_line(size =grid_line))+scale_y_continuous(limits=c(0,40),breaks=c(0,12,24,36))+geom_point(color=point_color,size=point_size)+stat_summary(fun.y=median,geom="point",size=1,aes(group=1,color='red'))
mtpt<-ggplot(a, aes(x=Clade2, y=pt.))+geom_violin()+coord_flip()+theme(axis.ticks.x = element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position = "none",panel.grid.minor.x = element_line(size =grid_line),panel.grid.minor.y = element_line(size =grid_line),panel.grid.major.x = element_line(size =grid_line),panel.grid.major.y = element_line(size =grid_line))+scale_y_continuous(limits=c(0,12),breaks=c(0,4,8,12))+geom_point(color=point_color,size=point_size)+stat_summary(fun.y=median,geom="point",size=1,aes(group=1,color='red'))
rep<-ggplot(a, aes(x=Clade2, y=repe.))+geom_violin()+coord_flip()+theme(axis.ticks.x = element_blank(),axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position = "none",panel.grid.minor.x = element_line(size =grid_line),panel.grid.minor.y = element_line(size =grid_line),panel.grid.major.x = element_line(size =grid_line),panel.grid.major.y = element_line(size =grid_line))+scale_y_continuous(limits=c(0,50),breaks=c(0,15,30,45))+geom_point(color=point_color,size=point_size)+stat_summary(fun.y=median,geom="point",size=1,aes(group=1,color='red'))
ggarrange( size, gc, genes, nrow=1, ncol=3 )
ggarrange( intron, rep, mtpt, nrow=1, ncol=3 )
dev.off()
