import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

#------------------dynamic data-----------------------------------------------------
# Load the dataset (replace with the path to your file)
file_path='/home/surface10/mrasolon/file_storage_git/CNP2_Best/MLacc_results.csv'
mlacc_results = pd.read_csv(file_path)
# storage of figures location
data_path='/home/surface10/mrasolon/file_storage_git/CNP2_Best/'
# mentionning the metric
metric=['R2','slope', 'dNRMSE'][2]
# mentionning the pool 
poolnames=['som', 'biomass','microbe', 'litter']
ipool=poolnames[3] #loop here

# mentionning the variable
name_prefix=['carbon', 'nitrogen', 'phosphorus']

var_labels=["Cpool", "Npool", "Ppool"]


#------------------static data-----------------------------------------------------

#ivar=name_prefix[1] #index
for index, ivar in enumerate(name_prefix):
   print(index)
   print(ivar) 
   # Filter for variable-related variables and group by 'var' and 'ind'
   meta_var_data=mlacc_results[mlacc_results['comp'].str.contains(ipool)]
   var_data = meta_var_data[meta_var_data['var'].str.contains(ivar)]

#-------------------------compartment dependences---------------------------------------------
   if ipool=='biomass':
     if index==0: 
          biomass_columns= ~meta_var_data['var'].str.endswith('_p') & ~meta_var_data['var'].str.endswith('_n') 
     if index==1: 
          biomass_columns=meta_var_data['var'].str.endswith('_n')
     if index==2:
          biomass_columns=meta_var_data['var'].str.endswith('_p')
      
     biomass_data=meta_var_data[biomass_columns]
     print(biomass_data)
     # Remove rows with missing values in 'var', 'ind', or 'R2'
     biomass_data_clean = biomass_data.dropna(subset=['var', 'ind', metric])

     # Create the pivot table
     heatmap_data = biomass_data_clean.pivot_table(index='ind', columns='var', values=metric).fillna(0)


     r2_2d = heatmap_data.values
     xTickLabel=["Leaf","SapAB","SapBE","HeartAB","HeartBE","Root","Fruit","Carbres","Labile"]
   
   if ipool=='litter':
     if index==0: 
          lit_columns= ~meta_var_data['var'].str.endswith('_p') & ~meta_var_data['var'].str.endswith('_n') 
     if index==1: 
          lit_columns=meta_var_data['var'].str.endswith('_n')
     if index==2:
          lit_columns=meta_var_data['var'].str.endswith('_p')

     lit_data=meta_var_data[lit_columns]
     heatmap_data = lit_data.pivot_table(index='var', columns='ind', values=metric) 
     temp = heatmap_data.values
     num_rows=temp.shape[0] // 2
     r2_2d=temp.reshape(num_rows, 6)
     xTickLabel=["StructuralAB","WoodyAB","MetaAB","StructuralBE", "WoodyBE", "MetaBE"]
   

   elif ipool=='microbe': 
      if index==0:
          mic_columns= ~meta_var_data['var'].str.endswith('_p') & ~meta_var_data['var'].str.endswith('_n')
      if index==1:
          mic_columns=meta_var_data['var'].str.endswith('_n')
      if index==2:
          mic_columns=meta_var_data['var'].str.endswith('_p')

      mic_data=meta_var_data[mic_columns]
      heatmap_data = mic_data.pivot_table(index='var', columns='ind', values=metric)
      temp = heatmap_data.values
      num_rows=temp.shape[0] // 2
      r2_2d=temp.reshape(num_rows, 4)
      xTickLabel=["RmicrAB","KmicrAB","RmicrBE","KmicrBE"]
   elif ipool=='som': 
     # Pivot the table so that 'var' becomes the row index and 'ind' becomes the column
      heatmap_data = var_data.pivot_table(index='var', columns='ind', values=metric)

      # Convert the pivoted table into a NumPy array
      r2_2d = heatmap_data.values
      xTickLabel=["Active","ChemProtect","PhysProtect"]
   #print(r2_2d)
#------------------------- end of compartment dependences------------------------------------
# Display the resulting array
#print(data_array)
 
   subLabel=var_labels[index]
 

   subps=len(xTickLabel)
   npfts=14
   yTickLabel = [
           "PFT02",
           "PFT03",
           "PFT04",
           "PFT05",
           "PFT06",
           "PFT07",
           "PFT08",
           "PFT09",
           "PFT10",
           "PFT11",
           "PFT12",
           "PFT13",
           "PFT14",
           "PFT15",
       ]
   yTickLabel = yTickLabel[0:npfts]
#----------------- colors of the maps setting-------------------------------------- 
   fonts = 7
   colors1 = plt.cm.YlGn(np.linspace(0, 1, 128))
   colors2 = plt.cm.YlGn_r(np.linspace(0, 1, 128))
   colors = np.vstack((colors1, colors2))
   mycolor_R2 = ["maroon", "tomato", "gold", "limegreen", "forestgreen"]
   mycolor_slope = [
         "maroon",
          "tomato",
         "gold",
           "limegreen",
           "forestgreen",
           "forestgreen",
           "limegreen",
           "gold",
           "tomato",
           "maroon",
       ]
   mycolor_rmse = ["forestgreen", "limegreen", "gold", "tomato", "maroon"]
   mymap = mcolors.LinearSegmentedColormap.from_list("my_colormap", colors)
   mymap_R2 = mcolors.LinearSegmentedColormap.from_list("my_list", mycolor_R2, N=5)
   mymap_slope = mcolors.LinearSegmentedColormap.from_list(        "my_list", mycolor_slope, N=10
        )
   mymap_rmse = mcolors.LinearSegmentedColormap.from_list("mylist", mycolor_rmse, N=5)
#-------------------------------- plotting ------------------------------------------
# matplotlib viz 
   fig, axs = plt.subplots(nrows=3, figsize=(8, 18))
   for jj in range(0, subps):
         # print(jj)
       for ii in range(0, npfts):
           # print(R22_n[ii,jj])
             # axs[0].text(-0.5 + jj, ii, str(r2_2d[ii, jj]), size=fonts, color="k")
           axs[0].text(-0.5 + jj, ii, f"{r2_2d[ii, jj]:.2f}", size=fonts, weight="bold", color="k")
   my_x_ticks = np.arange(subps)
   axs[0].set_xticks(my_x_ticks)
   # axs[0].set_xticklabels([""])
   axs[0].set_xticklabels(xTickLabel, rotation=60)
   my_y_ticks = np.arange(npfts)
   axs[0].set_yticks(my_y_ticks)
   axs[0].set_yticklabels(yTickLabel)
   axs[0].set_title(metric+"_" +subLabel) #index
   fig.subplots_adjust(right=0.9)
   l = 0.92
   b = 0.66
   w = 0.015
   h = 0.22
   rect = [l, b, w, h]
   cbar_ax = fig.add_axes(rect)
   sc = axs[0].imshow(r2_2d, vmin=0.5, vmax=1, cmap=mymap_R2)
   plt.colorbar(sc, cax=cbar_ax)
   plt.savefig(data_path+"Eval_"+metric+"_" + ipool + subLabel + ".png") #index
'''
#------------ the other metrics : slope and dNRMSE---------------
   # slope
   axs[1].imshow(slope_2d, vmin=0.75, vmax=1.25, cmap=mymap_slope)
   for jj in range(0, subps):
       for ii in range(0, npfts):
           axs[1].text(
               -0.5 + jj,
               ii,
               f"{slope_2d[ii, jj]:.2f}",
               size=fonts,
               color="k",
               weight="bold",
           )   
   my_x_ticks = np.arange(subps)
   axs[1].set_xticks(my_x_ticks)
   # axs[1].set_xticklabels([""])
   my_y_ticks = np.arange(npfts)
   axs[1].set_yticks(my_y_ticks)
   axs[1].set_yticklabels(yTickLabel)
   axs[1].set_title("slope_" + subLabel[0])
   fig.subplots_adjust(right=0.9)
   l = 0.92
   b = 0.39
   w = 0.015
   h = 0.22
   rect = [l, b, w, h]
   cbar_ax = fig.add_axes(rect)
   sc = axs[1].imshow(slope_2d, vmin=0.75, vmax=1.25, cmap=mymap_slope)
   plt.colorbar(sc, cax=cbar_ax)

# rmse
   axs[2].imshow(dNRMSE_2d, vmin=0, vmax=0.25, cmap=mymap_rmse)
   for jj in range(0, subps):
       for ii in range(0, npfts):
           axs[2].text(
               -0.5 + jj,
               ii,
               f"{dNRMSE_2d[ii, jj]:.2f}",
               size=fonts,
               color="k",
               weight="bold",
            )
   my_x_ticks = np.arange(subps)
   axs[2].set_xticks(my_x_ticks)
   axs[2].set_xticklabels(xTickLabel, rotation=60)
   my_y_ticks = np.arange(npfts)
   axs[2].set_yticks(my_y_ticks)
   axs[2].set_yticklabels(yTickLabel)
   axs[2].set_title("dNRMSE_" + subLabel[0])
   fig.subplots_adjust(right=0.9)
   l = 0.92
   b = 0.12
   w = 0.015
   h = 0.22
   rect = [l, b, w, h]
   cbar_ax = fig.add_axes(rect)
   sc = axs[2].imshow(dNRMSE_2d, vmin=0, vmax=0.25, cmap=mymap_rmse)
   plt.colorbar(sc, cax=cbar_ax)

   plt.savefig(data_path + "Eval_all_" + ipool + subLabel[0] + ".png")









#crreate the heatmap
plt.figure(figsize=(6, 8))
sns.heatmap(data_array, annot=True, cmap="RdYlGn", cbar_kws={'label': 'R²'}, 
            xticklabels=[1, 2, 3], yticklabels=range(2, 16))
plt.title("Eval_all_somCpool")
plt.xlabel("Index (ind)")
plt.ylabel("PFT")
output_filename = f"Eval_all"#_{component}{pool_label}.png"
plt.savefig(output_filename)
# Show the plot
plt.show()
'''
