import pandas as pd
import matplotlib.pyplot as plt

# directry
path = 'D:\\yasuda-lab\\B4\\save_img\\ODGMRF\\'

# csvs file name
file1 = path + 'denoised_for_DVGMRF' + '.csv'
file2 = path + 'origindata' + '.csv'

# read csv
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# x axis
x1 = df1['sample']

# y axis
y1 = df1['bit']
y2 = df2['bit']

# plot
plt.plot(x1, y1, label='DVGMRF')
plt.plot(x1, y2, color='red', label='original', alpha=0.675)

# plt.legend(fontsize='large')
plt.ylim(-1, 1)

plt.savefig(path + "\\" + 'compare' + ".png")
plt.show()