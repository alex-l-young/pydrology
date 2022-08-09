import pandas as pd
from datetime import datetime, timedelta

start_date = '2022-01-01'
end_date = '2022-01-10'

start_dt = datetime.strptime(start_date, '%Y-%m-%d')
end_dt = datetime.strptime(end_date, '%Y-%m-%d')

ts_a = []
ts_b = []
cur_dt = start_dt
while cur_dt <= end_dt:
    ts_a.append(cur_dt)
    ts_b.append(cur_dt)
    cur_dt = cur_dt + timedelta(days=1)

ts_a.pop(5)
ts_a.pop(6)
a = {
    'dog': [1, 2, 3, 4, 5, 7, 9, 10],
    'cat': [11, 12, 13, 14, 15, 17, 19, 20]
}

dfa = pd.DataFrame(a)
dfa['datetime'] = ts_a
print(dfa)

bseries = pd.Series(ts_b, name='datetime')

dfb = dfa.merge(bseries.to_frame(), on='datetime', how='right')
print(dfb)