#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 15:14:52 2018
@author: amayer
This script is sorting a list with the bubble method (values compared 2 by 2
and swapped if not in the right order)

"""

myList = [26, 54, 93, 17, 77, 31, 44, 55, 20]
N = len(myList)
i = 1

while (i < N-1):
    j = 0
    while (j < N-i):
        if myList[j] > myList[j+1]:
            temp = myList[j+1]
            myList[j+1] = myList[j]
            myList[j] = temp
            j += 1
        else:
            j += 1
    i += 1

print(myList)
