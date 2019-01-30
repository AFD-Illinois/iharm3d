#!/usr/bin/env python3

import sys
import pickle

avgs = []
for fname in sys.argv[1:]:
  avgs.append(pickle.load(open(fname, "rb")))

num_keys = [len(avg.keys()) for avg in avgs]
avg_w_max_keys = num_keys.index(max(num_keys))

uni = {}
for key in avgs[avg_w_max_keys].keys():
  uni[key] = avgs[avg_w_max_keys][key]
  for avg in avgs[1:]:
    if key in avg:
      # Keep track of averages w/weights, otherwise just sum since everything's time-dependent
      if key[-2:] == '_r' or key[-3:] == '_th':
        uni[key] += avg[key]*avg['avg_w']
      else:
        uni[key] += avg[key]

pickle.dump(uni, open("eht_out.p", "wb"))
