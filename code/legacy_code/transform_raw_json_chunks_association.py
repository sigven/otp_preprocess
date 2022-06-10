#!/usr/bin/env python

import os
import subprocess
letters = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1



i = 0
k = 1
while i < len(letters):
   j = 0
   #while j < 1:
   while j < len(letters):
      suffix = str(letters[i]) + str(letters[j])
      fname = 'target_validation_association_' + str(suffix)
      fname_pp = 'association_data_' + str(k) + '.json'
      if os.path.exists(fname):
         num_entries = file_len(fname)
         f = open(fname,'r')
         f2 = open(fname_pp,'w')
         f2.write("[\n")
         v = 1
         #print(str(num_entries))
         for line in f:
             #print(str(v))
             if v < num_entries:
                 #print(str(v) + '\t' + str(num_entries))
                 f2.write(line.rstrip() + ",\n")
             else:
                 f2.write(line.rstrip() + "\n")
             v = v + 1
         f2.write("]")
         f2.close()
         f.close()

         os.system('rm -f ' + str(fname))
         os.system('gzip ' + str(fname_pp))
         print(str(fname_pp))
         k = k + 1
      j = j + 1
   i = i + 1
