
import pandas as pd
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

import numpy as np
import seaborn as sns;

# print(plt.style.available)

def main():
    suffix = "titer_time.csv"
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/"
    output_dir = input_dir + "plots/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)



    data_RV = pd.read_csv(input_dir +"/RVB14/" + suffix)
    data_CV = pd.read_csv(input_dir +"/CVB3/" + suffix)
    data_all = pd.concat([data_RV, data_CV], sort=False)
    data_all["Virus"] = data_all["Virus"].apply(lambda x: x.split(" ")[0])
    # print(data_all.to_sting())

    data_all.to_csv(input_dir + "/titer_time_all.csv", sep=',', encoding='utf-8')


    plot1 = sns.catplot(x="Passage", y="PFU/mL", data=data_all, kind="point", hue="Virus")
    plot1.set(yscale='log')
    plot1.set(ylim=(10**6, 10**9))
    # plt.show()
    plt.savefig(output_dir + "titer_VS_time.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()