#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 10:16:27 2020

@author: daniel
"""


search_s=[]
m=0

for i in range(1000):
    print(i)
    m=1
    d=0
    while d!=9:
        
        m+=1
        rset,envus=get_mfs(reactions, ex_mets, phenotype, m)
        d = len(rset)
    search_s.append(m)