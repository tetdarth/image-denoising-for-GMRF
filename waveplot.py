import matplotlib.pyplot as plt
import os
import glob
import pandas as pd

path = "D:\\yasuda-lab\\B4\\save_img\\ODGMRF"
files = glob.glob(path + "\\*.csv")

for i in files:
    data = pd.read_csv(i)

    bit = data['bit'].astype('float64')
    title = data.columns[2]

    bit.plot()
    plt.title(title)
    plt.ylim(-1, 1)
    plt.savefig(path + "\\" + title + ".png")
    plt.show()
    os.remove(i)