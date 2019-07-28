import os, sys

import numpy as np
from astropy.table import Table

def read_catalog_sextractor(filename):
    f=open(filename)

    namelist=[]

    headflag=1
    while headflag:
        string=f.readline()
        if len(string)==0:
            break
        if not string[0]=='#':
            headflag=0
        else:
            name=string.split()[2]
            namelist.append(name)

    f.close()

    df=pd.read_csv(filename,delim_whitespace=1,comment='#',names=namelist)
    return df.copy()

def read_ps_data(filename):
    index_all = []
    values = []

    f_psdata = open(filename, 'r')

    line = f_psdata.readline() #First object header
    index = 0
    line = f_psdata.readline() #Column names
    names = line[:-1].split(',')
    line = f_psdata.readline() #Column format
    line = f_psdata.readline() #First row data

    if line[:13]=='no rows found':
        pass
    else:
        index_all.append(index)
        value = line[:-1].split(',')
        if len(value)==len(names):
            values.append(value.copy())

    while len(line)>0:
        line = f_psdata.readline()
        if line[:5]=='Input':
            index+=1
        elif line[:13]=='no rows found':
            continue
        else:
            value = line[:-1].split(',')
            if len(value)==len(names):
                index_all.append(index)
                values.append(value.copy())

    values = np.array(values)
    data = Table(values, names=names)

    for column in data.columns:
        try:
            data[column] = np.array(data[column], dtype=float)
        except ValueError:
            data[column] = np.array(data[column], dtype=str)

    data['OrigNum'] = index_all
    return data
