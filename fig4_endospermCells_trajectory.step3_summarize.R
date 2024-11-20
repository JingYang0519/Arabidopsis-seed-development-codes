library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(gplots)
library(viridis)

getwd()

work_path = "endo0514"
rds <-readRDS("refined_anno0514.endosperm.stage1-9.for_celltrajectory.rds")

time_point <- unique(rds$Organ)
time_point

replication_times=500
res_median_umap = list()
for(time_i in 1:(length(time_point)-1)){
  print(paste0(time_point[time_i], ":", time_point[time_i+1]))
  dat = readRDS(paste0(work_path, "/",time_point[time_i],"_",time_point[time_i+1],"_Knn_pca_permutation.rds"))
  state_1 = row.names(dat[[1]])
  state_1 = paste0(time_point[time_i+1], ":", gsub(paste0(time_point[time_i+1], ":"), "", state_1))
  state_2 = names(dat[[1]])
  state_2 = paste0(time_point[time_i], ":", gsub(paste0(time_point[time_i], ":"), "", state_2))
  tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))
  for(i in 1:nrow(dat[[1]])){
    for(j in 1:ncol(dat[[1]])){
      xx = NULL
      for(k in 1:replication_times){
        xx = c(xx, dat[[k]][i,j])
      }
      tmp_1[i,j] = median(xx[!is.na(xx)])
    }
  }
  tmp_1 = data.frame(tmp_1)
  row.names(tmp_1) = state_1
  names(tmp_1) = state_2
  res_median_umap[[time_i]] = tmp_1
}

dat = NULL
for(i in 1:length(res_median_umap)){
  print(time_point[i])
  dat = rbind(dat, melt(as.matrix(res_median_umap[[i]])))
}

dat = data.frame(dat)
names(dat) = c("nex", "pre", "prob")

saveRDS(dat, paste0(work_path, "/edge_all.rds"))

### Here we use “cell state” to mean an annotated cluster at a given stage. 

print(paste0("how many edges: ", nrow(dat)))
print(paste0("how many edges (> 0): ", nrow(dat[dat$prob>0,])))
print(paste0("how many edges (> 0.5): ", nrow(dat[dat$prob>=0.5,])))
print(paste0("how many edges (> 0.7): ", nrow(dat[dat$prob>=0.7,])))
print(paste0("how many edges (> 0.8): ", nrow(dat[dat$prob>=0.8,])))
print(paste0("how many nodes: ", length(unique(c(as.vector(dat$pre), as.vector(dat$nex))))))
print(paste0("how many cell types: ", length(unique(c(as.vector(dat$pre_cell), as.vector(dat$nex_cell))))))


#### extract edges with prob > 0.4
x = dat[dat$prob>=0,]
x = x[,c("pre","nex","prob")]
print(paste0("how many nodes now: ", length(unique(c(as.vector(x$pre), as.vector(x$nex))))))

res = x


dat_sub = res
dat_sub$pre_cell = unlist(lapply(as.vector(dat_sub$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat_sub$nex_cell = unlist(lapply(as.vector(dat_sub$nex), function(x) strsplit(x,"[:]")[[1]][2]))

sum(dat_sub$pre_cell == dat_sub$nex_cell)
sum(dat_sub$pre_cell == dat_sub$nex_cell)/nrow(dat_sub)
sum(dat_sub$pre_cell != dat_sub$nex_cell)
sum(dat_sub$pre_cell != dat_sub$nex_cell)/nrow(dat_sub)


print(paste0("how many edges: ", nrow(res)))

print(paste0("how many nodes: ", length(unique(c(as.vector(res$pre), as.vector(res$nex))))))

nex_list = as.vector(unique(res$nex))
tree = NULL

for(i in 1:length(nex_list)){
  res_sub = res[res$nex==nex_list[i],]

  if(nrow(res_sub)==1){
    tree = rbind(tree, res_sub)
  } else {

    res_sub = res_sub[order(res_sub$prob, decreasing = TRUE),]
    tree = rbind(tree, res_sub[1,])
  }
}

tree = data.frame(tree)
tree = tree[,c("pre","nex")]

write.table(res, paste0(work_path, "/edge_prob.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(tree, paste0(work_path, "/endo_edge.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
