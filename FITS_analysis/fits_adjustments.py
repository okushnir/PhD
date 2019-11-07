import os
import numpy as np
import pandas as pd


def fits_adjustments(data):
    data.sort_values(by=['pos'])
    grouped = data.groupby(['pos', 'gen'])
    df_all = pd.DataFrame(columns=["gen", "base", "freq", "pos"])
    for group in grouped:
        df = pd.concat([group[1]]*2, ignore_index=True)
        df.at[0, 'freq'] = 1 - df.at[1, 'freq']
        df.at[0, 'base'] = 0
        df.at[1, 'base'] = 1
        df_all = pd.concat([df_all, df], ignore_index=True)
    df_all["pos"] = df_all["pos"].astype(int)
    return df_all

def main():
    data = pd.read_csv(
        "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/190627_RV_CV/merged/RVB14/fits/test/p0-p12/SYN_MU~1.CSV")
    data = fits_adjustments(data)
    data.sort_values(by=['pos'])
    print(data.to_string())

if __name__ == "__main__":
    main()
# print(grouped.first())


# print(data.to_string())