# Script 1/2 to ensure that corr did the same thing between R and Python (June 2016)

#library(print)

METHOD<-'spearman'
my_mat<-matrix(c(3,5,1,2,4,2,9,NaN,8,1,NaN,3,5,2,8),nrow=5,ncol=3)
cat(sprintf("matrix:\n"))
prmatrix(my_mat)

# Compute correlation for matrix
out_cor<-cor(my_mat,use='pairwise.complete.obs',method=METHOD)
cat(sprintf("out_cor:\n"))
prmatrix(out_cor)

# Expect col1 v. col2 to be [3,5,2,4] v. [2,9,8,1]
col1_v_col2_cor<-cor(c(3,5,2,4),c(2,9,8,1),method=METHOD)
cat(sprintf("1 v. 2: %s, %s\n", col1_v_col2_cor, all.equal(col1_v_col2_cor,out_cor[1,2])))

# Expect col1 v. col3 to be [5,1,2,4] v. [3,5,2,8]
col1_v_col3_cor<-cor(c(5,1,2,4),c(3,5,2,8),method=METHOD)
cat(sprintf("1 v. 3: %s, %s\n", col1_v_col3_cor, all.equal(col1_v_col3_cor,out_cor[1,3])))

# Expect col2 v. col3 to be [9,8,1] v. [3,2,8]
col2_v_col3_cor<-cor(c(9,8,1),c(3,2,8),method=METHOD)
cat(sprintf("2 v. 3: %s, %s\n", col2_v_col3_cor, all.equal(col2_v_col3_cor,out_cor[2,3])))
