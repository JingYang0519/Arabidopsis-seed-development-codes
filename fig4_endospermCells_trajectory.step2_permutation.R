library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(htmlwidgets)
library(plotly)
library(monocle3)
library(FNN)
source("Fig5_endospermCells_trajectory.help_code.R")

work_path = "endo0517"
rds <-readRDS("refined_anno0514.endosperm.stage1-9.for_celltrajectory.rds")
time_point <- unique(rds$Organ)
time_point

for(kk in 1:(length(time_point)-1)){
    #kk = as.numeric(1)
    #rds <-readRDS("refined_anno0514.endosperm.stage1-9.for_celltrajectory.rds")
    time_i = time_point[kk]
    time_j = time_point[kk+1]
    print(time_i)
    print(time_j)
    emb = readRDS(paste0(work_path, "/", time_i, "_", time_j, "_umap3.rds"))
    emb = data.frame(emb)
    print(dim(emb))

    if(ncol(emb) != 3){
        print(XXX)
    }

   
    anno1 = subset(rds, Organ == time_i)
    Idents(anno1) <- "celltype0515"
    anno1$Anno = as.vector(anno1$celltype0515)
    state_counts1 <- table(anno1$Anno)
    states_to_remove1 <- names(state_counts1[state_counts1 < 3])
    print(time_i)
    print(states_to_remove1)
    cells_to_keep1 <- WhichCells(anno1, idents = setdiff(unique(Idents(anno1)), states_to_remove1))
    # 
    anno1 <- subset(anno1, cells = cells_to_keep1)
    # 
    anno1 <- AddMetaData(anno1, metadata = rep("pre", length(cells_to_keep1)), col.name = "day")
    anno1 <- AddMetaData(anno1, metadata = time_i, col.name = "stage")
    # 
    anno2 = subset(rds, Organ == time_j)
    Idents(anno2) <- "celltype0515"
    anno2$Anno = as.vector(anno2$celltype0515)
    state_counts2 <- table(anno2$Anno)
    states_to_remove2 <- names(state_counts2[state_counts2 < 3])
    print(time_j)
    print(states_to_remove2)
    cells_to_keep2 <- WhichCells(anno2, idents = setdiff(unique(Idents(anno2)), states_to_remove2))
    #
    anno2 <- subset(anno2, cells = cells_to_keep2)
    #
    anno2 <- AddMetaData(anno2, metadata = rep("pre", length(cells_to_keep2)), col.name = "day")
    anno2 <- AddMetaData(anno2, metadata = time_j, col.name = "stage")

    # 
    merge.data = merge(anno1, anno2)

    # 
    anno1 = subset(merge.data , Organ == time_i)
    anno1$Anno = as.vector(anno1$celltype0515)
    anno1@meta.data$day <- "pre"
    anno1 = anno1[[]][,c("day", "Anno")]
    anno1$stage = time_i

    anno2 = subset(merge.data , Organ == time_j)
    anno2$Anno = as.vector(anno2$celltype0515)
    anno2@meta.data$day <- "nex"
    anno2 = anno2[[]][,c("day", "Anno")]
    anno2$stage = time_j


    permutation_times = 1000
    k_neigh = 5

    res = list()
    for(rep_i in 1:permutation_times){
        
        anno1$state = anno1$Anno[sample(1:nrow(anno1))]
        anno2$state = anno2$Anno[sample(1:nrow(anno2))]
        
        anno = rbind(anno1, anno2)
        if(nrow(emb) != nrow(anno)){
            print("Error!")
            print(xxx)
        }
        pd = anno[rownames(emb),]
        
        emb_sub = emb
        pd_sub = pd
        
        irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
        irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
        pd_sub1 <- pd_sub[pd_sub$day == "pre",]
        pd_sub2 <- pd_sub[pd_sub$day == "nex",]
        #print(table(as.vector(pd_sub1$state)))
        pre_state_min = min(table(as.vector(pd_sub1$state)))
        paste("The pre_state_min is",pre_state_min)
        if (pre_state_min < k_neigh & pre_state_min >= 3){
            k_neigh = pre_state_min
            print(paste("The updated k_neigh is",k_neigh))
        }
        
        if (pre_state_min < 3){
            next
        }
        
        neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
        
        tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
        for(i in 1:k_neigh){
            tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
        }
        state1 <- names(table(as.vector(pd_sub1$state)))
        state2 <- names(table(as.vector(pd_sub2$state)))
        
        tmp2 <- matrix(NA,length(state2),length(state1))
        for(i in 1:length(state2)){
            x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
            for(j in 1:length(state1)){
                tmp2[i,j] <- sum(x==state1[j])
            }
        }
        tmp2 <- tmp2/apply(tmp2,1,sum)
        tmp2 <- data.frame(tmp2)
        row.names(tmp2) = state2
        names(tmp2) = state1
        
        res[[rep_i]] = tmp2
        
    }

    saveRDS(res, paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap_permutation.rds"))
}

x = list()
z = NULL
for(cnt in 1:(length(time_point)-1)){
    time_i = time_point[cnt]
    time_j = time_point[cnt+1]
    
    dat = readRDS(paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap_permutation.rds"))
    
    permutation_times = 1000
    y = NULL
    
    for(i in 1:permutation_times){
        y = c(y, as.vector(as.matrix(dat[[i]])))
    }
    
    x[[cnt]] = y
    z = c(z, y)
    
    print(paste0(time_i, ":", sum(y > 0.2)/length(y)))
}

print(sum(z >= 0.5)/length(z))

library(ggplot2)
dat = data.frame(edge_weights_by_permutation = z)

p<-ggplot(dat, aes(x=edge_weights_by_permutation)) +
    geom_histogram(position="identity", alpha=0.5, binwidth=0.01) + 
    geom_vline(xintercept = 0.5, colour = "red") +
    theme_classic(base_size = 15)
p

# 
ggsave(filename = "Knn_umap_edge_weights_by_permutation.pdf", plot = p, width = 8, height = 6, units = "in")

# 
ggsave(filename = "Knn_umap_edge_weights_by_permutation.png", plot = p, width = 8, height = 6, units = "in")
