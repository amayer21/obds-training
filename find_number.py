#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:03:39 2018
Find if the number 'x' is present in list of numbers 'L'
@author: amayer
"""


import sys


L = [26, 54, 93, 17, 77, 31, 44, 55, 20]
found = False
position = -1
index = 0
x = 18
x = int(sys.argv[1])

while found is False and (index <= len(L)-1):
    if L[index] == x:
        found = True
        position = index + 1
    else:
        index += 1
if found is True:
    print("The number {} is in the list at position {}".format(x, position))
else:
    print("The number {} isn't in the list".format(x))
