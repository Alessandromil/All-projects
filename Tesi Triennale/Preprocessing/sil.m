clc 
close all 

n = nan*ones(8,1);

%%solo HF, LZ3, STV 
X_30s = [[HF_sani.HF_30_33; n], [LZ3_sani.LZ3_30_33; n], [STV_sani.STV_30_33; n]]; 
X_30m = [HF_malati.HF_30_33, LZ3_malati.LZ3_30_33, STV_malati.STV_30_33]; 
X = [X_30s, X_30m]; 
clust = kmeans(X, 2); 

%30-33
silhouette(X,clust)
title('30-33')
close all 

X_34s = [[HF_sani.HF_34_36; n], [LZ3_sani.LZ3_34_36; n], [STV_sani.STV_34_36; n]]; 
X_34m = [HF_malati.HF_34_36, LZ3_malati.LZ3_34_36, STV_malati.STV_34_36]; 
X = [X_34s, X_34m]; 
clust = kmeans(X, 2); 

%34-36
silhouette(X,clust)
title('34-36')

X_37s = [[HF_sani.HF_37_38; n], [LZ3_sani.LZ3_37_38; n], [STV_sani.STV_37_38; n]]; 
X_37m = [HF_malati.HF_37_38, LZ3_malati.LZ3_37_38, STV_malati.STV_37_38]; 
X = [X_37s, X_37m]; 
clust = kmeans(X, 2); 

%37-38
silhouette(X,clust)
title('37-38')
