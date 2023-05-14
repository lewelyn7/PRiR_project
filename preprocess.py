import pandas as pd

data = pd.read_csv("results.csv", delimiter=";", index_col=None).drop(columns=['result'])

data = data.groupby(['mode', 'size', 'threads']).mean()

simple = data.loc[1]['time'] * 1000
parallel = data.loc[2]['time'] * 1000

data = pd.DataFrame({"simple": simple, "parallel": parallel})
data.to_csv("result_preprocessed.csv", decimal=',')


biggest_size = data.index[-1][0]
one_thread_value = data.loc[biggest_size].loc[2]['simple']

# for speedup chart, take the biggest size and print time for 1,2,4,8 threads
only_big_parallel = data.loc[biggest_size].drop(columns=['simple'])
only_big_parallel = pd.concat([only_big_parallel, pd.DataFrame([{"parallel": one_thread_value}], index=[1])])
only_big_parallel.sort_index(inplace=True)
only_big_parallel.to_csv("results_speedup.csv")

# for each size take the 1 threaded (simple) and parallel (8 threads) time
only_times = data.reorder_levels(['threads', 'size']).loc[8]
only_times.to_csv("results_only_max_threads.csv")


print("done :)")