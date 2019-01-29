#!/usr/bin/env python3

import sys
import pickle

avgs = []
for fname in sys.argv[1:]:
  avgs.append(pickle.load(open(fname, "rb")))

uni = {}
for key in avgs[0].keys():
  uni[key] = avgs[0][key]
  for avg in avgs[1:]:
    uni[key] += avg[key]

pickle.dump(uni, open("eht_out.p", "wb"))
