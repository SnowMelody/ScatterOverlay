import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import errstate, isneginf, array
import time
import argparse


# command-line interface
parser = argparse.ArgumentParser(
    description='ScatterOverlay. Generates overlapping scatter from 2 scatter plots.'
)
parser.add_argument('-size', type=float, default=0.0125, help='Marker size')
parser.add_argument('-colours', type=str, default='seagreen,rebeccapurple,orange', help='Scatter plot colour')
parser.add_argument('-samples', type=str, default='1,2,3,4', help='Samples to be plotted')
parser.add_argument('-filename', type=str, default='data.csv', help='CSV file to be read')
args = parser.parse_args()

theta = args.size
colourlist = [str(args.colours.split(',')[0]), str(args.colours.split(',')[1]), str(args.colours.split(',')[2])]
samplelist = [int(args.samples.split(',')[0]), int(args.samples.split(',')[1]), int(args.samples.split(',')[2]), int(args.samples.split(',')[3])]

print('Generating scatter plots...')
print('Plots will require some time to generate for files with large datasets.')
start_loop = time.time()

# import data and store in dataframe
data = pd.read_csv(args.filename)   
df = pd.DataFrame(data) 
df = df.set_index('Unnamed: 0')
df.index.name = None
 

# For first scatter
df_temp = pd.DataFrame(df.iloc[:,samplelist[0]-1])
df_temp = df_temp.join(df.iloc[:,samplelist[1]-1])
df_temp['Names'] = df_temp.index

# log10 transform (ignore 0)
with errstate(divide='ignore'):
    result = np.log10(df_temp.iloc[:,0] + 1)
    result2 = np.log10(df_temp.iloc[:,1] + 1)
result[isneginf(result)] = -1
result2[isneginf(result2)] = -1
df_temp.iloc[:,0] = result
df_temp.iloc[:,1] = result2

# remove col with neg values
df_temp = df_temp.loc[(df_temp.iloc[:,0]>=0)]
df_temp = df_temp.loc[(df_temp.iloc[:,1]>=0)]


# For second scatter
df_temp2 = pd.DataFrame(df.iloc[:,samplelist[2]-1])
df_temp2 = df_temp2.join(df.iloc[:,samplelist[3]-1])
df_temp2['Names'] = df_temp2.index

# log10 transform (ignore 0)
with errstate(divide='ignore'):
    result = np.log10(df_temp2.iloc[:,0] + 1)
    result2 = np.log10(df_temp2.iloc[:,1] + 1)
result[isneginf(result)] = -1
result2[isneginf(result2)] = -1
df_temp2.iloc[:,0] = result
df_temp2.iloc[:,1] = result2

# remove col with neg values
df_temp2 = df_temp2.loc[(df_temp2.iloc[:,0]>=0)]
df_temp2 = df_temp2.loc[(df_temp2.iloc[:,1]>=0)]


# find overlapping points, defined by a threshold value
x_val_inc = []  
y_val_inc = []
x_val_exc1 = []
y_val_exc1 = []
x_val_exc2 = []
y_val_exc2 = []
exc1_names = []
exc2_names = []

# combine both dataframes, then set threshold
for i in range(0, len(df_temp)):
    radius = np.sqrt((df_temp.iloc[i,0] - df_temp2.iloc[:,0])**2 + (df_temp.iloc[i,1] - df_temp2.iloc[:,1])**2)

    if (radius < (2*theta)).any() == True:
        x_val_inc.append(df_temp.iloc[i,0])
        y_val_inc.append(df_temp.iloc[i,1])
    else:
        x_val_exc1.append(df_temp.iloc[i,0])
        y_val_exc1.append(df_temp.iloc[i,1])
        exc1_names.append(df_temp.iloc[i,2])

for i in range(0, len(df_temp2)):
    radius = np.sqrt((df_temp2.iloc[i,0] - df_temp.iloc[:,0])**2 + (df_temp2.iloc[i,1] - df_temp.iloc[:,1])**2)

    if (radius < (2*theta)).any() == True:
        x_val_inc.append(df_temp2.iloc[i,0])
        y_val_inc.append(df_temp2.iloc[i,1])
    else:
        x_val_exc2.append(df_temp2.iloc[i,0])
        y_val_exc2.append(df_temp2.iloc[i,1])
        exc2_names.append(df_temp2.iloc[i,2])

    
df_inc = pd.DataFrame(x_val_inc)
df_inc.columns = ['X']
df_inc['Y'] = y_val_inc

df_exc = pd.DataFrame(x_val_exc1)
df_exc.columns = ['X']
df_exc['Y'] = y_val_exc1

df_exc2 = pd.DataFrame(x_val_exc2)
df_exc2.columns = ['X']
df_exc2['Y'] = y_val_exc2

# find names of differential genes
df_similar = pd.DataFrame(exc1_names)
df_similar2 = pd.DataFrame(exc2_names)

df_diff_genes = pd.DataFrame(np.concatenate((df_similar, df_similar2), axis=0))
df_diff_genes = df_diff_genes.drop_duplicates(subset=None, keep='first', inplace=False)

# write files
write_filename = args.filename.split('.csv')[0]
df_diff_genes.to_csv('%s_diff_genes.csv'%write_filename, index=False)
print('File written')

 
# plot actual figures
fig, ax = plt.subplots(1, 1, num='Scatter 1', figsize=(6,6))
ax.set_title(df_temp.columns[0] + '/' + df_temp.columns[1] + ' Scatter')
ax.set_xlabel(df_temp.columns[0])
ax.set_ylabel(df_temp.columns[1])
pts = plt.scatter(x=df_temp.iloc[:,0], y=df_temp.iloc[:,1], c=colourlist[0], marker='o', s=2)
plt.gca().set_aspect('equal')

fig2, ax2 = plt.subplots(1, 1, num='Scatter 2', figsize=(6,6))
ax2.set_title(df_temp2.columns[0] + '/' + df_temp2.columns[1] + ' Scatter')
ax2.set_xlabel(df_temp2.columns[0])
ax2.set_ylabel(df_temp2.columns[1])
pts2 = plt.scatter(x=df_temp2.iloc[:,0], y=df_temp2.iloc[:,1], c=colourlist[1], marker='o', s=2)
plt.gca().set_aspect('equal')

fig3, ax3 = plt.subplots(1, 1, num='Overlap', figsize=(6,6))
ax3.set_title('Overlapping Scatter')
ax3.set_xlabel('Log(X)')
ax3.set_ylabel('Log(Y)')
pts_inc = plt.scatter(x=df_inc.iloc[:,0], y=df_inc.iloc[:,1], c=colourlist[2], marker='o', s=2, zorder=1)
pts_exc = plt.scatter(x=df_exc.iloc[:,0], y=df_exc.iloc[:,1], c=colourlist[0], marker='o', s=2, zorder=2)
pts_exc2 = plt.scatter(x=df_exc2.iloc[:,0], y=df_exc2.iloc[:,1], c=colourlist[1], marker='o', s=2, zorder=2)
plt.gca().set_aspect('equal')

# make pie chart
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d} ({p:.2f}%)'.format(v=val, p=pct)
    return my_autopct

pc_values = [len(df_diff_genes), len(df_temp)-len(df_diff_genes)]
pc_labels = ['DE genes', 'Non-DE genes']

fig4, ax4 = plt.subplots(1, 1, num='Pie Chart', figsize=(6,6))
ax4.pie(x=pc_values, labels=pc_labels, autopct=make_autopct(pc_values), shadow=True)
ax4.axis('equal')  


# set axes limits so scale unchanged during dynamic plotting (also set all to be same for comparison)
x_min, x_max, y_min, y_max = ax.axis()
ax.set_xlim(-0.1, x_max)
ax.set_ylim(-0.1, y_max)
ax.xaxis.set_ticks(np.arange(0, 6, 1))
ax.yaxis.set_ticks(np.arange(0, 6, 1))

x2_min, x2_max, y2_min, y2_max = ax2.axis()
ax2.set_xlim(-0.1, x_max)
ax2.set_ylim(-0.1, y_max)
ax2.xaxis.set_ticks(np.arange(0, 6, 1))
ax2.yaxis.set_ticks(np.arange(0, 6, 1))

x3_min, x3_max, y3_min, y3_max = ax3.axis()
ax3.set_xlim(-0.1, x_max)
ax3.set_ylim(-0.1, y_max)
ax3.xaxis.set_ticks(np.arange(0, 6, 1))
ax3.yaxis.set_ticks(np.arange(0, 6, 1))


### Marker point converter ###
# set some variables
r = theta
N = 1

# Calculate radius in pixels (wrt each figure):
rr_pix = (ax.transData.transform(np.vstack([r, r]).T) -
          ax.transData.transform(np.vstack([np.zeros(N), np.zeros(N)]).T))
rpix, _ = rr_pix.T

rr_pix2 = (ax2.transData.transform(np.vstack([r, r]).T) -
          ax2.transData.transform(np.vstack([np.zeros(N), np.zeros(N)]).T))
rpix2, _ = rr_pix2.T

rr_pix3 = (ax3.transData.transform(np.vstack([r, r]).T) -
          ax3.transData.transform(np.vstack([np.zeros(N), np.zeros(N)]).T))
rpix3, _ = rr_pix3.T

# Calculate and update size in points (wrt each figure):
size_pt = (2*rpix/fig.dpi*72)**2
size_pt2 = (2*rpix2/fig2.dpi*72)**2
size_pt3 = (2*rpix3/fig3.dpi*72)**2

pts.set_sizes(size_pt)
pts2.set_sizes(size_pt2)
pts_inc.set_sizes(size_pt3)
pts_exc.set_sizes(size_pt3)
pts_exc2.set_sizes(size_pt3)


end_loop = time.time()
print('Scatter plots generated')
print('Total time taken:', end_loop - start_loop)


# save and display plots
fig.savefig('%s_scatter1.png'%write_filename)
fig2.savefig('%s_scatter2.png'%write_filename)
fig3.savefig('%s_overlay.png'%write_filename)
fig4.savefig('%s_piechart.png'%write_filename)
plt.show()
