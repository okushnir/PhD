
import pandas as pd
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

import numpy as np
import seaborn as sns;

# print(plt.style.available)

def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    output_dir = input_dir + "plots/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_titter = pd.read_csv(input_dir + "titter_time.csv")
    print(data_titter.to_string())
    plot1 = sns.catplot(x="Passage", y="PFU/mL", data=data_titter, kind="point")
    plt.savefig(output_dir + "titterVStime.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()