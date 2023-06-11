import pandas as pd
import sys
input_file_prefix = sys.argv[1]
data = pd.read_csv(input_file_prefix+".csv", delimiter=";", index_col=None).drop(columns=['result'])

data = data.groupby(['mode', 'size', 'threads']).mean()

simple = data.loc[1]['time'] * 1000
parallel = data.loc[2]['time'] * 1000

data = pd.DataFrame({"simple": simple, "parallel": parallel})
data.to_csv(input_file_prefix+"_preprocessed.csv", decimal=',')


biggest_size = data.index[-1][0]
one_thread_value = data.loc[biggest_size].loc[2]['simple']

# for speedup chart, take the biggest size and print time for 1,2,4,8 threads
only_big_parallel = data.loc[biggest_size].drop(columns=['simple'])
only_big_parallel = pd.concat([only_big_parallel, pd.DataFrame([{"parallel": one_thread_value}], index=[1])])
only_big_parallel.index.rename("threads", inplace=True)
only_big_parallel.sort_index(inplace=True)
only_big_parallel.to_csv(input_file_prefix+"_speedup.csv")

# for each size take the 1 threaded (simple) and parallel (8 threads) time
only_times = data.reorder_levels(['threads', 'size']).loc[8]
only_times.to_csv(input_file_prefix+"_only_max_threads.csv")


print("done :)")