from ctypes import string_at
import operator
from typing import final
from shapely.geometry import mapping, Polygon
from operator import itemgetter
import fiona
from geomet import wkt
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import geopandas as gpd
from geopandas import GeoDataFrame
import shapefile
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.geometry import MultiLineString
from shapely.geometry import Point
from road_lines_polygons import safety_toronto, safety_toronto_int, clustered_toronto, original_toronto_without, grid , original_toronto_with, accidents
import statistics
import math
import xlsxwriter

def areameanmedian(list_areas):
    s = 0
    for i in list_areas:
        s+=i
    mean = s/len(list_areas)
    median = statistics.median(list_areas)
    #print("this is mean: ", mean)
    #print("this is medeian: ", median)

    s = 0
    for l in list_areas:
        d = (l - mean) ** 2
        s += d
    paro = len(list_areas)
    std = round(math.sqrt(s/paro), 4)
    #print("this is std: ", std)
    return [mean, median, std]

def list_to_dgf(listt,epgs):
    d = {'geometry':[]}
    for pp in listt:
        d['geometry'].append(pp)

    gdf = gpd.GeoDataFrame(d, crs = epgs)
    return gdf

#takes line dataset (each line was mean,med,max,min width values), accident dataset -->
#associate for each line the accidents that is related with
# 1 lists --> length and accident number 
def acccheck(lns, pnts, epgs):
    acc = gpd.read_file(pnts)
    lns = gpd.read_file(lns)
    acc.crs = epgs
    lll = []
    lenght_accidents = []
    mnwidth_accidents = []
    medianwidth_accidents = []
    # For each line
    hn = []
    tr = []
    mr = []
    for a in range(len(lns)):
        # find accident (points) of each line
        rec = lns.loc[a,'geometry']
        rec2 = rec.buffer(9)
        n = acc.sindex.query(rec2, predicate='intersects')
        #If at least one accident occur in line
        if len(n) > 0:
            if rec.length <= 100:
                hn.append(1)
            elif rec.length > 100 and rec.length <= 300:
                tr.append(1)
            else:
                mr.append(1)
            lenght_accidents.append([rec.length, len(n)])
            mnwidth_accidents.append([lns.loc[a,'mean'], len(n)])
            if lns.loc[a,'median'] == 'only one line in the cluster':
                lns.loc[a, 'min'] =  np.float64(lns.loc[a, 'mean'])
                lns.loc[a, 'max'] =  np.float64(lns.loc[a, 'mean'])
                lns.loc[a,'median'] = np.float64(lns.loc[a,'mean'])
            else:
                lns.loc[a, 'min'] =  np.float64(lns.loc[a, 'min'])
                lns.loc[a, 'max'] =  np.float64(lns.loc[a, 'max'])
                lns.loc[a,'median'] = np.float64(lns.loc[a,'median'])
            medianwidth_accidents.append([lns.loc[a,'median'], len(n)])

            mn = lns.loc[a,'mean']
            md = lns.loc[a, 'median']
            minn = lns.loc[a, 'min']
            maxx = lns.loc[a, 'max']
            #print(minn)
            #print(maxx)
            num_accidents_permeter =  round(len(n) / rec.length, 3)
            
            lll.append([rec, len(n),rec.length, mn, md, minn, maxx, num_accidents_permeter] )

    print(len(hn))
    print(len(tr))
    print(len(mr))
    return lll, lenght_accidents

def writelineswithaccidents(lins):
    alllins = lins
    schema = {
        'geometry': 'LineString',
        'properties': {'accidents_number':'int', 'length_of_cluster':'float:13.3', 'mean':'float:13.3', 'median':'float:13.3','min':'float:13.3', 'max':'float:13.3', 'accident_norm': 'float:13.3'},
    }

    with fiona.open('lines_accidents.shp', 'w', 'ESRI Shapefile', schema) as c:
        ## If there are multiple geometries, put the "for" loop here
        for ii in alllins:
            c.write({
                'geometry': mapping(ii[0]),
                'properties': {'accidents_number':ii[1], 'length_of_cluster':ii[2], 'mean':ii[3], 'median':ii[4], 'min':ii[5], 'max':ii[6], 'accident_norm':ii[7]},
            })

def testt(pols2, lins2, epgs):
    print("START")
    lins = gpd.read_file(lins2)
    pols = gpd.read_file(pols2)
    l_final = []
    lins.crs = epgs
    mnwidth_normalized_accid = []
    mdwidth_normalized_accid = []
    for ii in range(len(lins)):
        rec = lins.loc[ii,'geometry']
        n = pols.sindex.query(rec, predicate = 'intersects')
        #print("line: ")
        #print(rec)
        if len(n) > 1:
           # print("more than 1 polygon")
            ind = 0
            for fsdf in n:
                ind+= pols.loc[fsdf, 'accident_p']
            final_ind = ind/len(n)
            #print(final_ind)
        elif len(n) == 1:
            #print(pols.loc[n[0], 'geometry'])
            final_ind = pols.loc[n[0], 'accident_p']
        else:
            print("no polygon intersect line")
            print(rec)
        
        l_final.append([lins.loc[ii,'geometry'], int(lins.loc[ii,'accidents_']), lins.loc[ii,'length_of_'], lins.loc[ii, 'mean'], lins.loc[ii, 'median'], lins.loc[ii, 'min'], lins.loc[ii, 'max'], lins.loc[ii, 'accident_n'], final_ind])
        
        rate_acc = lins.loc[ii, 'accident_n'] - final_ind
        mnwidth_normalized_accid.append([lins.loc[ii, 'mean'], rate_acc])
        mdwidth_normalized_accid.append([lins.loc[ii, 'median'],rate_acc])
        mn = lins.loc[ii,'mean']
        md = lins.loc[ii, 'median']
        minn = lins.loc[ii, 'min']
        maxx = lins.loc[ii, 'max']


    return l_final, mnwidth_normalized_accid, mdwidth_normalized_accid

def writefinallines(lins):
    alllins = lins
    schema = {
        'geometry': 'LineString',
        'properties': {'accidents_number':'float:13.3', 'length_of_cluster':'float:13.3', 'mean':'float:13.3', 'median':'float:13.3','min':'float:13.3', 'max':'float:13.3', 'accident_n': 'float:13.3', 'normalized_traf':'float:13.3'},
    }

    with fiona.open('lines_accidents_final.shp', 'w', 'ESRI Shapefile', schema) as c:
        ## If there are multiple geometries, put the "for" loop here
        for ii in alllins:
            c.write({
                'geometry': mapping(ii[0]),
                'properties': {'accidents_number':ii[1], 'length_of_cluster':ii[2], 'mean':ii[3], 'median':ii[4], 'min':ii[5], 'max':ii[6], 'accident_n':ii[7], 'normalized_traf':ii[8]},
            })

def stt(lins2):
    lins = gpd.read_file(lins2)
    mn = []
    md = []
    small_c = []
    l = []
    for i in range(len(lins)):
        rec_mn = lins.loc[i, 'mean']
        rec_md = lins.loc[i, 'median']
        rec_l = lins.loc[i, 'geometry'].length
        mn.append(rec_mn)
        if rec_md == 'only one line in the cluster':
            rec_md = np.float64(rec_mn)
        else:
            rec_md = np.float64(rec_md)
        md.append(rec_md)
        l.append(rec_l)
        if rec_l < 10.1:
            small_c.append(1)
    mm = areameanmedian(mn)[0]
    md = areameanmedian(md)[1]
    ml = areameanmedian(l)[0]
    mdl = areameanmedian(l)[1]
    sc = len(small_c)
    return [mm, md, ml, mdl, sc]


epgs = "EPSG:32617"

aaaa = stt(clustered_toronto)
print(aaaa)
new_lists = acccheck(original_toronto_with, accidents, epgs)
new_list = new_lists[1]

print("write excel for length and accidents")
with xlsxwriter.Workbook('test.xlsx') as workbook:
    worksheet = workbook.add_worksheet()

    for row_num, data in enumerate(new_list):
        worksheet.write_row(row_num, 0, data)

#write the final lines with the normalized accidents
lll = new_lists[0]
writelineswithaccidents(lll)


final_ls = testt(grid, "lines_accidents.shp", epgs)
#final_l = final_ls[0]
#writefinallines(final_l)


neo_list = final_ls[2]
with xlsxwriter.Workbook('neo_test2.xlsx') as workbook:
    worksheet = workbook.add_worksheet()

    for row_num, data in enumerate(neo_list):
        worksheet.write_row(row_num, 0, data)



# Not used !
road_polygons = safety_toronto 

def analysisroadsafety(lins2):
    lins = gpd.read_file(lins2)
    r1_mn = []
    r2_mn =[]
    r3_mn = []
    r4_mn = []
    r1_md = []
    r2_md =[]
    r3_md = []
    r4_md = []

    for ii in range(len(lins)):
        print("GO")
        if lins.loc[ii, 'mean'] <= 5:
            print(1)
            rate = lins.loc[ii,'accident_n'] - lins.loc[ii, 'normalized']
            r1_mn.append(rate)
        elif lins.loc[ii,'mean'] > 5 and lins.loc[ii,'mean'] <= 10:
            print(2)
            rate = lins.loc[ii,'accident_n'] - lins.loc[ii, 'normalized']
            r2_mn.append(rate)
        elif lins.loc[ii,'mean'] > 10 and lins.loc[ii,'mean'] <= 15:
            print(3)
            rate = lins.loc[ii,'accident_n'] - lins.loc[ii, 'normalized']
            r3_mn.append(rate)
        elif lins.loc[ii,'mean'] > 15:
            print(4)
            rate = lins.loc[ii,'accident_n'] - lins.loc[ii, 'normalized']
            r4_mn.append(rate)
        
        print(rate)

        if lins.loc[ii, 'median'] <= 5:
            print(1)
            rate1 = lins.loc[ii,'accident_n'] - lins.loc[ii, 'normalized']
            r1_md.append(rate1)
        elif lins.loc[ii,'median'] > 5 and lins.loc[ii,'median'] <= 10:
            print(2)
            rate1 = lins.loc[ii,'accident_n'] - lins.loc[ii, 'normalized']
            r2_md.append(rate1)
        elif lins.loc[ii,'median'] > 10 and lins.loc[ii,'median'] <= 15:
            print(3)
            rate1 = lins.loc[ii,'accident_n'] - lins.loc[ii, 'normalized']
            r3_md.append(rate1)
        elif lins.loc[ii,'median'] > 15:
            print(4)
            rate1 = lins.loc[ii,'accident_n'] - lins.loc[ii, 'normalized']
            r4_md.append(rate1)
        print(rate1)
    
    print(sum(r1_mn))
    print(sum(r2_mn))
    print(sum(r3_mn))
    print(sum(r4_mn))

    print(sum(r1_md))
    print(sum(r2_md))
    print(sum(r3_md))
    print(sum(r4_md))

def analysisroadsafety2(road_polygons, clust, lins, epgs):
    roads = gpd.read_file(road_polygons)
    acc_lins = gpd.read_file(lins)
    clust_lins = gpd.read_file(clust)
    acc_lins.crs = epgs
    clust_lins.crs = epgs
    pls = []
    for i in range(len(roads)):
        rec = roads.loc[i, 'geometry']
        n = acc_lins.sindex.query(rec, predicate='intersects')
        ratee = 0
        if len(n) > 0:
            print(rec)
            for ln in n:
                r = acc_lins.loc[ln,'accident_n']
                r1 = acc_lins.loc[ln, 'normalized']
                rat = r-r1
                print(rat)
                ratee += rat
        print(ratee)
        if ratee != 0:
            pls.append([rec, ratee])
    r1 =[]
    r2 =[]
    r3 = []
    for oo in pls:
        p = oo[0]
        print(p)
        print(oo[1])
        nn = clust_lins.sindex.query(p, predicate='intersects')
        if len(nn) <= 2:
            r1.append(oo[1])
        elif len(nn) > 2 and len(nn)<=5:
            r2.append(oo[1])
        elif len(nn) > 5:
            r3.append(oo[1])
    
    print(sum(r1))
    print(sum(r2))
    print(sum(r3))

def accidents(acc_shapefile, pol_shp, epgs):
    acc = gpd.read_file(acc_shapefile)
    pols = gpd.read_file(pol_shp)
    pols.crs = epgs
    big_lines_shp = gpd.read_file(big_lines_shp)
    original_lines_shp = gpd.read_file(original_lines_shp)
    polygons_with_accident = []
    for aa in range(len(acc)):
        rec = acc.loc[aa, 'geometry']
        print(rec)
        p_n1 = pols.sindex.query(rec, predicate="contains")
        print("length of polygons that contain the accident: ", len(p_n1))
        if len(p_n1) > 0:
            for pp in p_n1:
                plg = pols.loc[pp, 'geometry']
                plg_mean = pols.loc[pp, 'mean_w(m)']
                plg_type = pols.loc[pp, 'Road_part']
                plg_median = pols.loc[pp, 'median_w(m']
                print(plg)
                polygons_with_accident.append([plg, plg_mean, plg_median, plg_type])
        else:
            print("no polygon contains the accident, lets check for polygons that are close: ")
            rec_ex = rec.buffer(25)
            print(rec_ex)
            p_n2 = pols.sindex.query(rec_ex, predicate="intersects")
            print("length of polygons that are close to the accident: ", len(p_n2))
            if len(p_n2)>0:
                distances = []
                for ppp in range(len(p_n2)):
                    plg = pols.loc[ppp, 'geometry']
                    print(plg)
                    d = plg.boundary.distance(rec)
                    distances.append([plg,d, ppp])

                d_s = sorted(distances, key= operator.itemgetter(1))
                print("this is polygon with sortest distance to the accident: ", d_s[0][0])
                plg_mean = pols.loc[d_s[0][2], 'mean_w(m)']
                plg_type = pols.loc[d_s[0][2], 'Road_part']
                plg_median = pols.loc[d_s[0][2], 'median_w(m']
                polygons_with_accident.append([d_s[0][0], plg_mean, plg_median, plg_type])


            else:
                print("this point has no polygons close!")
                continue
        

        #b_n1 = big_lines_shp.sindex.query(rec, predicate="intersects")
        #o_n1 = original_lines_shp.sindex.query(rec, predicate="intersects")

def accidents2(acc_shapefile, pol_shp, epgs):
    acc = gpd.read_file(acc_shapefile)
    pols = gpd.read_file(pol_shp)
    pols.crs = epgs
    polygons_with_accident = []
    allpointsused = []
    if acc.crs != pols.crs:
        print("GAMW TIN TUXH MOU")
    for aa in range(len(pols)):
        rec = pols.loc[aa, 'geometry']
        rec_mean = pols.loc[aa, 'mean_w(m)']
        rec_type = pols.loc[aa, 'Road_part']
        rec_median = pols.loc[aa, 'median_w(m']
        rec_max = pols.loc[aa, 'max_w(m)']
        rec_min = pols.loc[aa, 'min_w(m)']
        rec_ex = rec.buffer(6)
        #print(rec_ex)
        p_n1 = acc.sindex.query(rec_ex, predicate="intersects")
        #print(len(p_n1))
        #print("length of accidents that the polygon is related with: ", len(p_n1))
        l = len(p_n1)
        if len(p_n1) > 0:
            for pp in p_n1:
                pnt = acc.loc[pp, 'geometry']
                print(pnt)
                if pnt not in allpointsused:    
                    allpointsused.append(pnt)
                else:
                    l -= 1
        else:
            continue
        
        if l>0:
            polygons_with_accident.append([rec, rec_mean, rec_median,rec_max, rec_min, rec_type, l])

    return polygons_with_accident

def accidentslines(acc_shapefile, original_lines_shp, big_lins, epgs):
    acc = gpd.read_file(acc_shapefile)
    orig_lins = gpd.read_file(original_lines_shp)
    big_lins = gpd.read_file(big_lins)
    orig_lins.crs = epgs
    big_lins.crs = epgs
    liness = []
    allpointsused = []
    if acc.crs != orig_lins.crs:
        print("GAMW TIN TUXH MOU")
    for ee in range(len(orig_lins)):
        rec = orig_lins.loc[ee, 'geometry']
        rec_mean = orig_lins.loc[ee, 'mean_w(m)']
        rec_median = orig_lins.loc[ee, 'median_w(m']
        rec_max = orig_lins.loc[ee, 'max_w(m)']
        rec_min = orig_lins.loc[ee, 'min_w(m)']
        rec_ex = rec.buffer(9)
        if rec_mean > 0.0:
            #print(rec_ex)
            p_n1 = acc.sindex.query(rec_ex, predicate="intersects")
            #print(len(p_n1))
            #print("length of accidents that the polygon is related with: ", len(p_n1))
            l = len(p_n1)
            if len(p_n1) > 0:
                for pp in p_n1:
                    pnt = acc.loc[pp, 'geometry']
                    #print(pnt)
                    if pnt not in allpointsused:    
                        allpointsused.append(pnt)
                    else:
                        l -= 1
            else:
                continue
            
            if l>0:
                liness.append([rec, rec_mean, rec_median, rec_max, rec_min, l])
        else:
            continue
    
    linesss = []
    allpointsused2 = []
    for a in range(len(big_lins)):
        rec = big_lins.loc[a, 'geometry']
        rec_mean = big_lins.loc[a, 'mean_w(m)']
        rec_median = big_lins.loc[a, 'median_w(m']
        rec_max = big_lins.loc[a, 'max_w(m)']
        rec_min = big_lins.loc[a, 'min_w(m)']
        rec_ex = rec.buffer(9)
        if rec_mean > 0.0:
            #print(rec_ex)
            p_n1 = acc.sindex.query(rec_ex, predicate="intersects")
            #print(len(p_n1))
            #print("length of accidents that the polygon is related with: ", len(p_n1))
            l = len(p_n1)
            if len(p_n1) > 0:
                for pp in p_n1:
                    pnt = acc.loc[pp, 'geometry']
                    #print(pnt)
                    if pnt not in allpointsused2:    
                        allpointsused2.append(pnt)
                    else:
                        l -= 1
            else:
                continue
            
            if l>0:
                linesss.append([rec, rec_mean, rec_median, rec_max, rec_min, l])
        else:
            continue
    
    return [liness, linesss]


def writepolsacc(polygons_accidents):
    allpols = polygons_accidents
    schema = {
        'geometry': 'Polygon',
        'properties': {'mean_w(m)':'float:13.3', 'median_w(m)':'float:13.3', 'max_w(m)':'float:13.3', 'min_w(m)':'float:13.3', 'Road_part':'str', 'Accidents': 'int'},
    }

    with fiona.open('accident_polygons.shp', 'w', 'ESRI Shapefile', schema) as c:
        for i in allpols:
            c.write({
                'geometry': mapping(i[0]),
                'properties': {'mean_w(m)': i[1], 'median_w(m)': i[2], 'max_w(m)': i[3], 'min_w(m)' : i[4], 'Road_part' : i[5],'Accidents': i[6]},
            })

def writeoriginallinesacc(lines_accidents):
    lll = lines_accidents
    schema = {
        'geometry': 'LineString',
        'properties': {'mean_w(m)':'float:13.3', 'median_w(m)':'float:13.3', 'max_w(m)':'float:13.3', 'min_w(m)':'float:13.3', 'Accidents': 'int'},
    }

    with fiona.open('accident_original_lines.shp', 'w', 'ESRI Shapefile', schema) as c:
        for i in lll:
            c.write({
                'geometry': mapping(i[0]),
                'properties': {'mean_w(m)': i[1], 'median_w(m)': i[2], 'max_w(m)': i[3], 'min_w(m)' : i[4],'Accidents': i[5]},
            })

def writedissovledlinesacc(lines_acc_dis):
    lll = lines_acc_dis
    schema = {
        'geometry': 'LineString',
        'properties': {'mean_w(m)':'float:13.3', 'median_w(m)':'float:13.3', 'max_w(m)':'float:13.3', 'min_w(m)':'float:13.3', 'Accidents': 'int'},
    }

    with fiona.open('accident_dissolved_lines.shp', 'w', 'ESRI Shapefile', schema) as c:
        for i in lll:
            c.write({
                'geometry': mapping(i[0]),
                'properties': {'mean_w(m)': i[1], 'median_w(m)': i[2], 'max_w(m)': i[3], 'min_w(m)' : i[4],'Accidents': i[5]},
            })
