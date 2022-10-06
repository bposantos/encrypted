######################## import all libraries ###########################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D
import os
import argparse
import sys

######################## Args Parser ####################################
parser = argparse.ArgumentParser(description='3D plot/PCA/kmeans clustering')
parser.add_argument("-f", '--file', help="file_name.csv", action="store")
parser.add_argument("-m", '--method', help="the plot method: 'pca' or 'manual'", action="store")
parser.add_argument("-cl", '--n_clusters', help="file name", action="store")
parser.add_argument("-c1", '--column1', help="the columns positions to be added in the vizualization", action="store")
parser.add_argument("-c2", '--column2', help="the columns positions to be added in the vizualization", action="store")
parser.add_argument("-c3", '--column3', help="the columns positions to be added in the vizualization", action="store")
parser.add_argument("-o", '--output', help="file name", action="store")
args = parser.parse_args()

######################## Manipulate Dataframe ####################################

df2 = pd.read_csv(args.file) #backup file for csv generation

df1 = pd.read_csv(args.file) #file that will be manipulated
df1.drop(['Unnamed: 0','Amidated_Mass','[M+H+]','Peptides','Name'], axis=1, inplace=True)
cols = df1.columns.tolist()
cols = cols[-1:] + cols[:-1]
df1 = df1[cols]

#methods for data evaluation
#method1 - Unsupervised learning - PCA
#method2 - Physicochemical properties comparison

############################## Colors ########################################
'''
plot_color_gradients('Perceptually Uniform Sequential',
                     ['viridis', 'plasma', 'inferno', 'magma', 'cividis'])

plot_color_gradients('Sequential',
                     ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'])

plot_color_gradients('Sequential (2)',
                     ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                      'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
                      'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'])

plot_color_gradients('Diverging',
                     ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'])

plot_color_gradients('Cyclic', ['twilight', 'twilight_shifted', 'hsv'])

plot_color_gradients('Qualitative',
                     ['Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
                      'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b',
                      'tab20c'])

plot_color_gradients('Miscellaneous',
                     ['flag', 'prism', 'ocean', 'gist_earth', 'terrain',
                      'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap',
                      'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet',
                      'turbo', 'nipy_spectral', 'gist_ncar'])

'''
############################## PCA ###########################################

if args.method == 'pca':

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

   df2['Clusters'] = previsoes

   #Adding kmeans as target in the table
   df2.to_csv(args.output+'.csv')

   # import relevant libraries for 3d graph
   fig = plt.figure()
    
   # choose projection 3d for creating a 3d graph
   #axis = fig.add_subplot(1,1,1, projection='3d')
   #axis.set_xlabel("PCA 1", fontsize=12, weight='bold')
   #axis.set_ylabel("PCA 2", fontsize=12, weight='bold')
   #axis.set_zlabel("PCA 3", fontsize=12, weight='bold')

   u_labels = np.unique(previsoes)

   # x[:,0]is pc1,x[:,1] is pc2 while x[:,2] is pc3

   color_list = ['flag', 'prism','ocean', 'gist_earth', 'terrain',
   'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap',
   'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']

   #color_list =['twilight', 'twilight_shifted', 'hsv']

   print(len(color_list))

   for i in range(0,len(color_list)):
      axis = fig.add_subplot(4,4,i+1, projection='3d')
      scatter = axis.scatter(df[:,0],df[:,1],df[:,2],
         s=40,
         c=df2['Clusters'], 
         cmap=color_list[i],
         linestyles='solid',
         alpha=0.8, 
         linewidth=0.8, 
         edgecolor='black',
         picker=True)
      axis.set_xlabel("PCA 1", fontsize=4)
      axis.set_ylabel("PCA 2", fontsize=4)
      axis.set_zlabel("PCA 3", fontsize=4)

      #legend1 = axis.legend(*scatter.legend_elements(),
      #                    loc="upper left", title="Clusters")
      #axis.add_artist(legend1)

   plt.show()

   '''
   #Pick points
   def onpick(event):
      ind = event.ind[0]+1
      print(ind)

   fig.canvas.mpl_connect('pick_event', onpick)
   '''

   #show graphics
   #plt.show()

   ############################### Matrix of PCA components contribution ##############################################
   '''
   loadings = principal.components_
   #print('loadings',loadings)

   num_pc = principal.n_features_
   #print('num_pc',num_pc)

   pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
   #print('pc_list',pc_list)

   loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))

   #df2.drop(['Clusters'],axis=1,inplace=True)
   loadings_df['variable'] = df1.columns.values
   loadings_df = loadings_df.set_index('variable')
   #print(loadings_df)

   # get correlation matrix plot for loadings
   import seaborn as sns
   import matplotlib.pyplot as plt
   ax = sns.heatmap(loadings_df, annot=True, cmap='Spectral')
   plt.show()
   '''
############################## Manual columns selection ###########################################

if args.method == 'manual':            
   X = df1.loc[:,[args.column1,args.column2,args.column3]]

   scaling=StandardScaler()
    
   # Use fit and transform method
   scaling.fit(X)
   df=scaling.transform(X)

   from sklearn.cluster import KMeans
   kmeans = KMeans(n_clusters=int(args.n_clusters),random_state=0)
   previsoes = kmeans.fit_predict(df)

   df2['Clusters'] = previsoes

   #Adding kmeans as target in the table
   df2.to_csv(args.output+'.csv')

   fig=plt.figure()
   ax = fig.add_subplot(111, projection='3d')

   u_labels = np.unique(previsoes)
   
   scatter = ax.scatter(X.iloc[:,0],X.iloc[:,1],X.iloc[:,2],
      s=100,
      c=df1['Clusters'],
      cmap='cool',
      marker='o',
      linestyles='solid',
      alpha=0.8, 
      linewidth=1, 
      edgecolor='black',
      picker=True)
   
   ax.set_xlabel(args.column1, fontsize=10)
   ax.set_ylabel(args.column2, fontsize=10)
   ax.set_zlabel(args.column3, fontsize=10)

   legend1 = ax.legend(*scatter.legend_elements(),
                       loc="upper left", title="Clusters")
   ax.add_artist(legend1)

   #Pick points
   def onpick(event):
      ind = event.ind[0]+1
      print(ind)

   fig.canvas.mpl_connect('pick_event', onpick)

   plt.show()