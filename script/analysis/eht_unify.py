#!/usr/bin/env python3

import sys
import pickle
import numpy as np

avgs = []
for fname in sys.argv[1:]:
  print("Loading {}".format(fname))
  avgs.append(pickle.load(open(fname, "rb")))
  avgs[-1]['fname'] = fname

#for avg in avgs:
#  print("Name: {}, contents: {}".format(avg['fname'], avg.keys()))

num_keys = [len(avg.keys()) for avg in avgs]
avg_max_keys = num_keys.index(max(num_keys))

direct_list = ['fname', 'a', 'r', 'th', 'th_bh', 'th_5', 'th_bz', 'avg_start', 'avg_end']
keys_to_sum = [key for key in avgs[avg_max_keys].keys() if key not in direct_list]

uni = {}
for key in keys_to_sum:
  uni[key] = np.zeros_like(avgs[avg_max_keys][key])
  for avg in avgs:
    if key in avg:
      # Keep track of averages w/weights, otherwise just sum since everything's time-dependent
      if key[-2:] == '_r' or key[-3:] == '_th' or key[-4:] == '_rth' or key[-6:] == '_thphi':
        uni[key] += avg[key]*avg['avg_w']
      else:
        uni[key] += avg[key]

for key in direct_list:
  if key in avgs[avg_max_keys].keys():
    uni[key] = avgs[avg_max_keys][key]

with open("eht_out.p", "wb") as outf:
  print("Writing eht_out.p")
  pickle.dump(uni, outf)