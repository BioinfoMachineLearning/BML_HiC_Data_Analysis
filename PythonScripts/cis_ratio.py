import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import datetime

def plot_data(data, columns):

    df = pd.DataFrame(data=data, columns=columns)
    plt.figure(figsize=(10, 6))
    sns.barplot(x="section", hue="Sample", y="% of total cis", data=df)
    plt.title("Bar plot showing the % of cis reads(section/total cis)")
    plt.savefig("bars.png")

def generate():
    files = ["stats/27950_stats_dedup.pairsam", "stats/27951_stats_dedup.pairsam", "stats/27989_stats_dedup.pairsam", "stats/28033_stats_dedup.pairsam"]
    i = files[3]
    with open(i) as f:
        stor = {}
        for l in f:
            tmp = l.strip().split("\t")
            stor[tmp[0]] = int(tmp[1])

    long = stor["cis_20kb+"]
    cis  = stor["cis"]

    print("file is : "+ i)
    # >20kb
    long_percent = long/cis *100
    print("long: {:,} percentage is : ".format(long) + str(round(long_percent, 2)))


    # 1kb -20kb
    mid = stor["cis_1kb+"] - long
    mid_percent = mid/cis *100
    print("mid: {:,}  percentage is : ".format(mid) + str(round(mid_percent, 2)))

    #print(stor.keys())
    # <1kb
    close = cis - stor["cis_1kb+"]
    close_percent = close/cis *100
    print("close : {:,}  percentage is : ".format(close) + str(round(close_percent, 2)))

    print("Cis : {:,}".format(cis))

    print("Trans : {:,}".format(stor["trans"]))

columns=["section", "Sample", "% of total cis"]
data = [["close(<1kb)", "27950", 69.94], ["close(<1kb)", "27951", 71.31], ["close(<1kb)", "27989", 62.2], ["close(<1kb)", "28033", 59.42],
        ["mid(<1kb-20kb)", "27950", 7.21], ["mid(<1kb-20kb)", "27951", 6.6], ["mid(<1kb-20kb)", "27989", 8.77], ["mid(<1kb-20kb)", "28033", 10.28],
        ["far(>20kb)", "27950", 22.85], ["far(>20kb)", "27951", 22.09], ["far(>20kb)", "27989", 29.04], ["far(>20kb)", "28033", 30.3]]
plot_data(data=data, columns=columns)