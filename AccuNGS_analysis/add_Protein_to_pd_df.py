import pandas as pd
import numpy as np

def add_Protein_to_pd_df_func(pddf, region_lst):
    """

    :param pddf: pandas DataFrame of pipeline post manipulation
    :param region_lst: list of all the position number from the begging till the end. exaple [629, 835, 1621, 2329, 3196
                                                                                            , 3634, 3925, 4915, 5170,
                                                                                            5239, 5785, 7165]
    :return:
    """
    pddf["Protein"] = np.where(pddf["Pos"] <= region_lst[-12], "5'UTR",
                                             np.where(pddf["Pos"] <= region_lst[-11], "VP4",
                                                      np.where(pddf["Pos"] <= region_lst[-10], "VP2",
                                                               np.where(pddf["Pos"] <= region_lst[-9], "VP3",
                                                                        np.where(pddf["Pos"] <= region_lst[-8],
                                                                                 "VP1",
                                                                                 np.where(
                                                                                     pddf["Pos"] <= region_lst[-7],
                                                                                     "2A",
                                                                                     np.where(pddf[
                                                                                                  "Pos"] <= region_lst[-6], "2B",
                                                                                              np.where(
                                                                                                  pddf[
                                                                                                      "Pos"] <= region_lst[-5],
                                                                                                  "2C",
                                                                                                  np.where(
                                                                                                      pddf[
                                                                                                          "Pos"] <= region_lst[-4],
                                                                                                      "3A",
                                                                                                      np.where(
                                                                                                          pddf[
                                                                                                              "Pos"] <= region_lst[-3],
                                                                                                          "3B",
                                                                                                          np.where(
                                                                                                              pddf[
                                                                                                                  "Pos"] <= region_lst[-2],
                                                                                                              "3C",
                                                                                                              np.where(
                                                                                                                  pddf[
                                                                                                                      "Pos"] <= region_lst[-1],
                                                                                                                  "3D",
                                                                                                                  "3'UTR"))))))))))))
    return pddf