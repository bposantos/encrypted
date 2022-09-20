# import all libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
import argparse
import sys

######################## Args Parser ####################################
parser = argparse.ArgumentParser(description='3D plot/PCA/kmeans clustering')
parser.add_argument("-f", '--file', help="file_name.csv", action="store")
parser.add_argument("-o", '--output', help="file name", action="store")
parser.add_argument("-cl", '--n_clusters', help="file name", action="store")
args = parser.parse_args()

######################## Manipulate Dataframe ####################################

df1 = pd.read_csv(args.file)
del df1["Unnamed: 0"]
df1.drop(['Amidated_Mass','[M+H+]'], axis=1, inplace=True)
cols = df1.columns.tolist()
cols = cols[-1:] + cols[:-1]
df1 = df1[cols]
df1.set_index('Peptides', inplace=True) 

############################## PCA ###########################################

# Scale data before applying PCA
scaling=StandardScaler()
 
# Use fit and transform method
scaling.fit(df1)
Scaled_data=scaling.transform(df1)
 
# Set the n_components=3
principal=PCA(n_components=3)
principal.fit(Scaled_data)
df=principal.transform(Scaled_data)

print("The PC variance are: PC1 =", "%.2f" % principal.explained_variance_ratio_[0],\
   "PC2 =", "%.2f" % principal.explained_variance_ratio_[1],"PC3 =", "%.2f" % principal.explained_variance_ratio_[2])

print("The cumulative PCA", "%.2f" % principal.explained_variance_ratio_.cumsum()[2])

# Check the dimensions of data after PCA
print("The Matrix shape: ",df.shape)

################################## Clustering using kmeans #####################################################

from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=int(args.n_clusters),random_state=0)
previsoes = kmeans.fit_predict(df)

df1['Clusters']=previsoes

#Adding kmeans as target in the table
df1.to_csv(args.output+'.csv')

# import relevant libraries for 3d graph
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(10,10))
 
# choose projection 3d for creating a 3d graph
axis = fig.add_subplot(111, projection='3d')

u_labels = np.unique(previsoes)

# x[:,0]is pc1,x[:,1] is pc2 while x[:,2] is pc3
scatter = axis.scatter(df[:,0],df[:,1],df[:,2],
   s=100,
   c=df1['Clusters'], 
   cmap='gist_rainbow',
   linestyles='solid',
   alpha=0.8, 
   linewidth=1, 
   edgecolor='black', 
   picker=True)
axis.set_xlabel("PCA 1", fontsize=12, weight='bold')
axis.set_ylabel("PCA 2", fontsize=12, weight='bold')
axis.set_zlabel("PCA 3", fontsize=12, weight='bold')

legend1 = axis.legend(*scatter.legend_elements(),
                    loc="lower left", title="Clusters")
axis.add_artist(legend1)

#Pick points
def onpick(event):
   ind = event.ind[0]+1
   print(ind)

fig.canvas.mpl_connect('pick_event', onpick)

#show graphics
plt.show()

############################### Matrix of PCA components contribution ##############################################

loadings = principal.components_
#print('loadings',loadings)

num_pc = principal.n_features_
#print('num_pc',num_pc)

pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
#print('pc_list',pc_list)

loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
#print(loadings_df)

df1.drop(['Clusters'],axis=1,inplace=True)
loadings_df['variable'] = df1.columns.values
loadings_df = loadings_df.set_index('variable')
#print(loadings_df)

# get correlation matrix plot for loadings
import seaborn as sns
import matplotlib.pyplot as plt
ax = sns.heatmap(loadings_df, annot=True, cmap='Spectral')
plt.show()
