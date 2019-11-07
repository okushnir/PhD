import pandas as pd
import numpy as np
import glob
from Context_analysis_RV import checkKey


def main():
    virus = "CVB3"

    input_dir = "/Users/odedkushnir/PycharmProjects/Context/Runs/%s" % (virus)
    dirs = glob.glob(input_dir + "/%s*" % virus)
    pos_dic = {"RVB14": 3606, "CVB3": 3745}
    frame_pos = checkKey(pos_dic, virus)
    for passage in dirs:
        sample = glob.glob(passage + "/*.freqs")
        sample = sample[0]

    # sample = "/Users/odedkushnir/PycharmProjects/Context/RVB14-p2/RVB14-p2.freqs"

        data = pd.read_table(sample)
        data = data[data['Pos'] == np.round(data['Pos'])]
        data = data[data['Pos'] >= frame_pos]
        if virus == "CVB3":
            data = data[data['Pos'] <= 7299]
        data["Pos"] = data['Pos'].apply(lambda x: x-(frame_pos-1))

        data = data[data["Base"] != "-"]
        # data = data.drop(columns="Unnamed: 0")
        data.to_csv(sample, sep="\t", encoding='utf-8', index=False)

        with open(sample, 'r') as file:
            table = file.read()
            with open(sample, 'w') as file1:
                file1.write("#oded\n" + table)
                file1.close()


if __name__ == "__main__":
    main()
