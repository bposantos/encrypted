# import all libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#open table file
df1 = pd.read_csv('')
del df1["Unnamed: 0"]

# Scale data before applying PCA
scaling=StandardScaler()
 
# Use fit and transform method
scaling.fit(df1)
Scaled_data=scaling.transform(df1)
 
# Set the n_components=3
principal=PCA(n_components=3)
principal.fit(Scaled_data)
x=principal.transform(Scaled_data)

print("The PC variance are: PC1 =", "%.2f" % principal.explained_variance_ratio_[0],\
   "PC2 =", "%.2f" % principal.explained_variance_ratio_[1],"PC3 =", "%.2f" % principal.explained_variance_ratio_[2])

print("The cumulative PCA", "%.2f" % principal.explained_variance_ratio_.cumsum()[2])

# Check the dimensions of data after PCA
print("The Matrix shape: ",x.shape)

####################################################################################################################

# Matrix of PCA components contribution

loadings = principal.components_
#print(loadings)

num_pc = principal.n_features_
#print(num_pc)

pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
#print(pc_list)

loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
#print(loadings_df)

loadings_df['variable'] = df1.columns.values
loadings_df = loadings_df.set_index('variable')

print(loadings_df)

# get correlation matrix plot for loadings
import seaborn as sns
import matplotlib.pyplot as plt
ax = sns.heatmap(loadings_df, annot=True, cmap='Spectral')
plt.show()

####################################################################################################################

#Clustering using kmeans
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=3,random_state=0)
previsoes = kmeans.fit_predict(x)
df1['target']=previsoes

#Adding kmeans as target in the table
df1.to_csv('leish_df_kmeans.csv')

# import relevant libraries for 3d graph
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(10,10))
 
# choose projection 3d for creating a 3d graph
axis = fig.add_subplot(111, projection='3d')
 
# x[:,0]is pc1,x[:,1] is pc2 while x[:,2] is pc3
axis.scatter(x[:,0],x[:,1],x[:,2], c=df1['target'], cmap='plasma', picker=True)
axis.set_xlabel("PC1", fontsize=10)
axis.set_ylabel("PC2", fontsize=10)
axis.set_zlabel("PC3", fontsize=10)

#Pick points
def onpick(event):
   ind = event.ind[0]
   print(ind)

fig.canvas.mpl_connect('pick_event', onpick)

#show graphics
plt.show()
