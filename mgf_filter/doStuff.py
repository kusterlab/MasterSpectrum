from pyteomics import mgf
import csv
import pandas as pd
import numpy as np
from operator import itemgetter

class points(object):
    summed = 0
    count = 0

class Extractor(object):
    amzintensity = 0 

def extractSpectraFast(x):
    rows = 0
    with mgf.read(x) as spectra:
        for spectrum in spectra:
            rows = rows + len(spectrum['m/z array'])
    # first argument is rownumber, sec argument is columnnumber
    n = np.empty((rows,2))
    with mgf.read(x) as spectra:
        i = 0
        for spectrum in spectra:
            m = spectrum['m/z array']
            intensity = spectrum['intensity array']
            for m,intensity in zip(spectrum['m/z array'],spectrum['intensity array']):
                n[i,0] = m
                n[i,1] = intensity
                i = i + 1
    return np 

def extractSpectraFastForR(x):
    i = 1
    with mgf.read(x) as spectra, open("/home/tobiass/df.csv","wt") as csvfile:
        writr = csv.writer(csvfile)
        for spectrum in spectra:
            for m,intensity in zip(spectrum['m/z array'],spectrum['intensity array']):
                writr.writerow((m, intensity,i))
            i = i + 1

def calculateRelativeIntensity(aIntensity):
    relIntensity = np.array(aIntensity)
    relIntensity = relIntensity/np.amax(relIntensity)
    return relIntensity.tolist()

def extractSpectraWeighted(x):
    #df = pd.DataFrame({'mz':[],'intensity':[],'amount':[]})
    df = pd.DataFrame()
    count = 1
    loop = 0
    with mgf.read(x) as spectra:
        for spectrum in spectra:
            #print("##########")
            #print(count)
            count = count + 1
            df_2_add = pd.DataFrame({'mz':spectrum['m/z array'],
                                     'intensity':calculateRelativeIntensity(spectrum['intensity array']),
                                     'amount':1})
            df = df.append(df_2_add)
            for intensity_ordered in df['intensity'].order().drop_duplicates()[::-1]:
                #print("---------------------------------")
                ### exist ordered intensity 
                if ( (df['intensity'] == intensity_ordered ).any() ):
                    #print(intensity_ordered)
                    # problem: peaks with same intensity, which pick first
                    for int in df.query('intensity ==@intensity_ordered')['mz'].tolist():
                        #print(int)
                        border = calculateResolutionWindow(int)
                        #print(border)
                        expr = '@border[0] < mz < @border[1]'
                        pos2add = df.query(expr, parser='pandas')
                        loop = loop + 1
                        #print("loop unwichtig",loop)
                        add = pos2add.sum()
                        #print(add)
                        positionToDelete =  pos2add['mz'].tolist()
                        #print('length before: ',len(df.index))
                        df = df.query('mz not in @positionToDelete')
                        #print('length: ',len(df.index))
                        #### should i use summed rel. intensity or average rel.intensity
                        # scalar values must be frames to inserted into a Dataframe
                        df = df.append(pd.DataFrame({'mz'        : [add['mz']/add['amount']], 
                                                    'intensity' : [add['intensity']/add['amount']] ,
                                                     'amount': [add['amount']] }) , ignore_index=True)
            #print(len(df.index))
    #print(loop)
    return df

def sortRestrainByAmount(a,b,sortbyPosition = 0,restrainAmount = 0):
    sortIt = [list(x) for x in zip(*sorted(zip(a,b), key = itemgetter(sortbyPosition), reverse = True))]
    l = len(sortIt[0])
    if (l >restrainAmount):
        return [sortIt[0][0:restrainAmount],sortIt[1][0:restrainAmount]]
    else:
        return sortIt
    #filterIt = [(x,y) for x,y in zip(sortIt[0],sortIt[1]) if x > restrainvalue):
    

def extractSpectraWeightedByMathias(x):
    df = pd.DataFrame({'sum_weighted_mz':[],'sum_weighted_intensity':[],'count':[],'left_border':[],'right_border':[], 'startmz':[]})
    up = 1
    append_case = 0
    easy_merge = 0
    multiple_merge = 0
    with mgf.read(x) as spectra:
        for spectrum in spectra:
            print(up)
            up = up +1 
            mz_unsorted = spectrum['m/z array']
            intensity_unsorted = spectrum['intensity array']
            sortIt = sortRestrainByAmount(mz_unsorted,intensity_unsorted,1,100)
            for mz,intensity in zip(sortIt[0],sortIt[1]):
                search_in_area = df.query('left_border < @mz < right_border')
                entries_in_window = len(search_in_area.index)
                if(entries_in_window == 0):
                    border = calculateResolutionWindow(mz)
                    df = df.append(pd.DataFrame({'sum_weighted_mz' : [intensity * mz],
                                             'sum_weighted_intensity': [intensity],
                                             'left_border':[border[0]],
                                             'right_border': [border[1]],
                                             'count':[1],
                                             'startmz':mz}), ignore_index=True)
                    append_case = append_case + 1

                elif(entries_in_window == 1):
                    sum_rel_int = intensity + search_in_area['sum_weighted_intensity'].tolist()[0]
                    counter = search_in_area['count'].tolist()[0] + 1
                    sum_mz = search_in_area['sum_weighted_mz'].tolist()[0] + intensity * mz
                    position = sum_mz/sum_rel_int
                    border = calculateResolutionWindow(position)
                    pos_2_chg = search_in_area['sum_weighted_mz'].tolist()[0]
                    startm = search_in_area['startmz'].tolist()[0]
                    if (False):
                        print(mz, intensity)
                        print(search_in_area)
                        print(pos_2_chg)
                        print(sum_mz)
                        print(sum_rel_int)
                        print(counter)
                        print(border[0])
                        print(border[1])
                        print(df.loc[df.sum_weighted_mz == pos_2_chg])

                    #df.loc[df.sum_weighted_mz == pos_2_chg,] = [1,1,1,1,1]
                    df.loc[df.sum_weighted_mz == pos_2_chg,:] = [sum_mz, sum_rel_int, counter, border[0], border[1] ,startm  ]
                    easy_merge = easy_merge + 1

                else:
                    # old stuff
                    sum_rel_mz = sum(search_in_area['sum_weighted_mz']*search_in_area['sum_weighted_intensity'])
                    sum_rel_int = sum(search_in_area['sum_weighted_intensity'])
                    sum_counter = sum(search_in_area['count'])

                    sum_rel_mz = sum_rel_mz + mz * intensity
                    sum_rel_int = sum_rel_int + intensity
                    sum_counter = sum_counter + 1

                    position = sum_rel_mz/sum_rel_int
                    border = calculateResolutionWindow(position)
                    
                    #delete first
                    df = df.drop(search_in_area.index)

                    # now add
                    df = df.append(pd.DataFrame({'sum_weighted_mz' : [sum_rel_mz],
                                             'sum_weighted_intensity': [sum_rel_int],
                                             'left_border':[border[0]],
                                             'right_border': [border[1]],
                                             'count':sum_counter,
                                             'startmz':-1}), ignore_index=True)
                     
                    #print(mz, intensity)
                    #print(search_in_area)
                    print("is bloed")
                    multiple_merge = multiple_merge +1
                    #exit(1)
            print(append_case, easy_merge, multiple_merge)
            print("####################################")

    return df

def extractSpectraStatistics(x):
    df = pd.DataFrame({'sum_weighted_mz':[],'sum_weighted_intensity':[],'count':[],'left_border':[],'right_border':[], 'startmz':[]})
    up = 1
    with mgf.read(x) as spectra, open("/home/tobiass/data/count_peaks_start_mgf/df.csv","wt") as csvfile:
        writr = csv.writer(csvfile)
        for spectrum in spectra:
            rel_int = calculateRelativeIntensity(spectrum['intensity array'])
            rel_int.sort(reverse = True)
            df = pd.DataFrame({'intensity':rel_int })
            l = len(df.index)
            percent_1 = len(df[df.intensity>0.01])
            percent_01 = len(df[df.intensity>0.001])
            if (l >= 50):
                top_50 = df.iloc[49]['intensity']
                top_50_exist = 1
            else:
                top_50 = df.iloc[-1]['intensity']
                top_50_exist = 0 
            if (l >= 100):
                top_100 = df.iloc[99]['intensity']
                top_100_exist = 1
            else:
                top_100 = df.iloc[-1]['intensity']
                top_100_exist = 0 

            writr.writerow((l, top_100, top_100_exist, top_50, top_50_exist, percent_1,percent_01))
            print(up)
            up = up + 1

def calculateResolutionWindow(x):
    return [x-x/calculateResolution(x),x+x/calculateResolution(x)]

def calculateResolution(mz):
    return np.power(10,5.847 + np.multiply(np.log10(mz),-0.546))
