import pandas as pd

data = pd.read_csv("results.csv", delimiter=";", index_col=None).drop(columns=['result'])

data = data.groupby(['mode', 'size']).mean()

simple = data.loc[1]['time'] * 1000
parallel = data.loc[2]['time'] * 1000

data = pd.DataFrame({"simple": simple, "parallel": parallel})
data.to_csv("result_preprocessed.csv", decimal=',')
print(data.head())