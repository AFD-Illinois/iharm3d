#!/usr/bin/env python3

avgs = []
for fname in sys.argv[1:]:
  avgs.append(pickle.load(open(fname, "rb")))

for key in avgs[0].keys():
  uni[key] = avgs[0][key]
  for avg in avgs[1:]:
    uni[key] += avgs[key]

pickle.dump(uni, open("eht_out.p", "wb"))