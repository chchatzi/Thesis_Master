from re import M
import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from shapely.coords import CoordinateSequence
#from matplotlib import figure
from shapely.geometry import mapping, Polygon
import fiona
from geomet import wkt
import pandas as pd
import geopandas as gpd
from geopandas import GeoDataFrame
import shapefile
from shapely.geometry import box
from shapely.geometry.base import EMPTY
from shapely.validation import explain_validity
import networkx as nx
import osmnx as ox
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.geometry import MultiLineString
from shapely.geometry import Point
from road_lines_polygons import  case1, edgescase1,case_narrow, edges_narrow, edges_parking, case1_a,case1_b,case1_c,case1_d
from road_lines_polygons import   case2, case2_a, case2_b, case2_c, case2_d
from operator import itemgetter
from shapely.ops import linemerge
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
import os
from shapely.ops import linemerge
import statistics

edges = gpd.read_file(edges_parking)
ground_truth_case1 = gpd.read_file("C:/Users/chcha/Desktop/GEOMATICS/2nd Year/Thesis/Clustering_datasets/ground_truth_5cases/case1/ground_truth_final.shp")
ground_truth_case2 = gpd.read_file("C:/Users/chcha/Desktop/GEOMATICS/2nd Year/Thesis/Clustering_datasets/ground_truth_5cases/parking_case/ground_truth.shp")
pls_case = gpd.read_file(case2_c)
clust_approach = gpd.read_file("C:/Users/chcha/Desktop/GEOMATICS/2nd Year/Thesis/p2_technical/code/methodology1/clust.shp")


#Help functions
def areameanmedian(list_areas):
    s = 0
    for i in list_areas:
        s+=i
    mean = s/len(list_areas)
    median = statistics.median(list_areas)
    s = 0
    for l in list_areas:
        d = (l - mean) ** 2
        s += d
    paro = len(list_areas)-1
    std = round(math.sqrt(s/paro), 4)
    return [mean, median, std]

def width(prepedicularline, polygon, roadline, id):
    #prependicular line between 2 offsets intersects polygon
    ln = prepedicularline.intersection(polygon)
    #if returns 2 or more linestrings we need to find the one that intersects the roadline
    if ln.type == 'MultiLineString':
        ll = list(ln)
        s = 0
        for i in ll:
            if i.intersects(roadline):
                measureline = i
                tom = measureline.length
                w = round(tom, 4)
                ln = measureline
                #print(measureline)
            else:
                s+=1
        if len(ll) > 1 and s == len(ll):
            d = []
            dd = []
            print("WEIRD CASE")
            for i in ll:
                midi = i.interpolate(0.5, normalized = True)
                midroadline = roadline.interpolate(0.5,normalized = True)
                distancetoroadline = LineString([midi, midroadline]).length
                d.append([distancetoroadline,i])
                dd.append(distancetoroadline)
            mindist = min(dd)
            for val in d:
                if val[0] == mindist:
                    tom = val[1].length
                    w = round(tom,6)
                    ln = val[1]
                    #print(ln)
    else:
        tom = ln.length
        w = round(tom, 4)
        #print(ln)
    
    return [w, ln, id]

def clusters_num(pls, ground_truth, clust_approach):
    print(len(pls))
    roads_samenum_clust = []
    for i in range(len(pls)):
        rec = pls.loc[i, 'geometry']
        #print(rec)
        n = ground_truth.sindex.query(rec, 'intersects')
        n2 = clust_approach.sindex.query(rec, 'intersects')
        #print(len(n))
        if len(n) == len(n2):
            roads_samenum_clust.append(rec)
    print(len(roads_samenum_clust))
    print(len(roads_samenum_clust)*100/len(pls))
    #for pp in roads_samenum_clust:
        #print(pp)
    return (len(roads_samenum_clust)*100/len(pls))

def indicators(pls,ground_truth,clust_approach):
    polygon_gtlines = []
    polygon_clustlines = []
    #for each polygon, find the lines(clusters) of ground truth and of clustering approach
    #indicator1_2 = []
    #indicator1_1 = []
    #indicator1_0 = []
    for i in range(len(pls)):
        truth = []
        clust = []
        rec = pls.loc[i, 'geometry']
        print("For polygon: ", rec)
        n = ground_truth.sindex.query(rec, 'intersects')
        n2 = clust_approach.sindex.query(rec, 'intersects')
        for ind in n:
            ln = ground_truth.loc[ind, 'geometry']
            truth.append(ln)
        for ind2 in n2:
            ln2 = clust_approach.loc[ind2, 'geometry']
            clust.append(ln2)  
        #Indicator 1 --> Number of clusters
        if len(truth) == len(clust):
            print("same num of clusters: ", len(truth))
            indicator1 = 1
            print(indicator1)
        else:
            print("diff num of clusters: ", len(truth), "and", len(clust))
            if len(truth) < len(clust):
                indicator1 = round(len(truth)/len(clust),3)
                print(indicator1)
            else:

                indicator1 = round(len(clust)/len(truth),3)
                print(indicator1)
        #for indicator 2 and 3
        #make a list with lines (clusters) of ground truth and basic stats (length, mean/median withd) for each line(cluster)
        lines_ground_truth_with_stats = []
        idd = 0
        for ll in truth:
            #print(ll)
            lengg = round(ll.length, 2)
            how_many_measure = round(lengg)
            interval = 1
            measure_lines_width = []
            for hh in range(how_many_measure):
                m = LineString([ll.interpolate(interval-1, normalized = False), ll.interpolate(interval, normalized = False)])
                #print(m)
                interval+=1
                mdpoint = m.interpolate(0.5, normalized = True)
                leftline = m.parallel_offset(50, 'left')
                rightline = m.parallel_offset(50, 'right')
                midpointleft = leftline.interpolate(0.5, normalized = True)
                midpointright = rightline.interpolate(0.5, normalized = True)
                prependicular_line = LineString([midpointleft, midpointright])
                prep_line = [prependicular_line, 1]
                #print(prep_line[0])
                #for each prep line find the 2 closest points to the polygon that contains the initial line (create the measuring line)
                inn = pls.sindex.query(mdpoint, predicate="within")
                if len(inn)!= 1 and len(inn)!= 0:
                    print("error")
                if len(inn)==0:
                    continue
                pol_contain_line = pls.loc[inn[0], 'geometry']
                oww = width(prep_line[0], pol_contain_line, ll, prep_line[1])
                measure_lines_width.append(oww[0])
            s = areameanmedian(measure_lines_width)
            lines_ground_truth_with_stats.append([ll,lengg, round(s[0],3), round(s[1],3), idd])
            idd+=1
        #print("Ground truth results to: ", len(truth), " lines, with lengths, mean width and median width: ")
        polygon_gtlines.append([rec, len(truth), lines_ground_truth_with_stats, indicator1])
        
        #make a list with lines (clusters) of clustering approach and basic stats (length, mean/median withd) for each line(cluster)
        lines_clustering_intersectionwithoriginallines_with_stats = []
        idd2 = 0
        for ll2 in clust:
            #print(ll2)
            #print(ll)
            lengg2 = round(ll2.length, 2)
            how_many_measure2 = round(lengg2)
            interval2 = 1
            measure_lines_width2 = []
            for hh2 in range(how_many_measure2):
                m2 = LineString([ll2.interpolate(interval2-1, normalized = False), ll2.interpolate(interval2, normalized = False)])
                #print(m)
                interval2+=1
                mdpoint2 = m2.interpolate(0.5, normalized = True)
                leftline2 = m2.parallel_offset(50, 'left')
                rightline2 = m2.parallel_offset(50, 'right')
                midpointleft2 = leftline2.interpolate(0.5, normalized = True)
                midpointright2 = rightline2.interpolate(0.5, normalized = True)
                prependicular_line2 = LineString([midpointleft2, midpointright2])
                prep_line2 = [prependicular_line2, 1]
                #print(prep_line[0])
                #for each prep line find the 2 closest points to the polygon that contains the initial line (create the measuring line)
                inn2 = pls.sindex.query(mdpoint2, predicate="within")
                if len(inn2)!= 1 and len(inn2)!= 0:
                    print("error")
                if len(inn2)==0:
                    continue

                pol_contain_line2 = pls.loc[inn2[0], 'geometry']
                oww2 = width(prep_line2[0], pol_contain_line2, ll2, prep_line2[1])
                measure_lines_width2.append(oww2[0])
            s2 = areameanmedian(measure_lines_width2)
            lines_clustering_intersectionwithoriginallines_with_stats.append([ll2, lengg2, round(s2[0],3), round(s2[1],3), idd2])
            idd2+=1
        polygon_clustlines.append([rec, len(clust), lines_clustering_intersectionwithoriginallines_with_stats,indicator1])
    return  polygon_gtlines, polygon_clustlines

def indicators_comparison(polygon_gtlines, polygon_clustlines):
    for pol in polygon_clustlines:
        for cc in pol[2]:
            if cc[1] < 6:
                print("i removed:")
                print(cc)
                pol[2].remove(cc)
    for pol2 in polygon_gtlines:
        for cc2 in pol2[2]:
            if cc2[1] < 6:
                print("i removed:")
                print(cc2)
                pol2[2].remove(cc2)
    polygons_indicators = []
    for p in polygon_gtlines:
        for pp in polygon_clustlines:
            if p[0] == pp[0]:
                print(p[0])
                print("clustering approach has this num of clusters: ")
                print(len(pp[2]))
                #for each cluster of clustering approach find the cluster of ground truth that we should compare with (the one with the biggest intersection)
                #small_diff_l = []
                #small_diff_w = []
                for cls in pp[2]:
                    print
                    print("this is cluster of approach1")
                    print(cls[0])
                    lengths_intersections = []
                    inds2 = []
                    inds3 = []
                    for cls_truth in p[2]:
                        intersectionn = cls[0].intersection(cls_truth[0].buffer(1))
                        if intersectionn.geom_type == 'LineString' or intersectionn.geom_type == 'MultilineString' and (intersectionn.is_empty != True):
                            #print(intersectionn)
                            l = intersectionn.length
                            lengths_intersections.append([l,cls_truth[0], cls_truth[1], cls_truth[2], cls_truth[3]])
                        else:
                            continue
                    lis = sorted(lengths_intersections, key=lambda x: x[0])
                    correct_cluster_ground_truth = [lis[-1][1], lis[-1][2], lis[-1][3], lis[-1][4]]
                    print("this is the corresponding cluster of ground truth")
                    print(correct_cluster_ground_truth[0])
                    print(cls)
                    print(correct_cluster_ground_truth)
                    #Compare length for indicator2
                    if cls[1] < correct_cluster_ground_truth[1]:
                        indicator2_c = cls[1] / correct_cluster_ground_truth[1]
                        inds2.append(indicator2_c)
                    else:
                        indicator2_c = correct_cluster_ground_truth[1]/ cls[1]
                        inds2.append(indicator2_c)

                    print("this is he final index for this cluster for indicator 2:", indicator2_c)

                    #width difference depends of length
                    if cls[1] <= 110:
                        if cls[2] < correct_cluster_ground_truth[2]:
                            ind3_mn = round(cls[2] / correct_cluster_ground_truth[2], 3) - 0.1
                        else:
                            ind3_mn = round(correct_cluster_ground_truth[2] / cls[2], 3)- 0.1
                        print(ind3_mn)
                        if cls[3] < correct_cluster_ground_truth[3]:
                            ind3_md = round(cls[3]/correct_cluster_ground_truth[3], 3)- 0.1
                        else:
                            ind3_md = round(correct_cluster_ground_truth[3]/cls[3], 3)- 0.1
                        print(ind3_md)
                        inds3_final = (ind3_mn + ind3_md) / 2
                        print('this is the final index for this cluster for indicator 3: ', inds3_final)
                        inds3.append(inds3_final)

                    elif cls[1] > 110 and cls[1] <= 210:
                        if cls[2] < correct_cluster_ground_truth[2]:
                            ind3_mn = round(cls[2] / correct_cluster_ground_truth[2], 3) - 0.05
                            if ind3_mn > 1:
                                ind3_mn = 1
                        else:
                            ind3_mn = round(correct_cluster_ground_truth[2] / cls[2], 3) - 0.05
                            if ind3_mn > 1:
                                ind3_mn = 1
                        print(ind3_mn)
                        if cls[3] < correct_cluster_ground_truth[3]:
                            ind3_md = round(cls[3]/correct_cluster_ground_truth[3], 3) - 0.05
                            if ind3_md > 1:
                                ind3_md = 1
                        else:
                            ind3_md = round(correct_cluster_ground_truth[3]/cls[3], 3) - 0.05
                            if ind3_md > 1:
                                ind3_md = 1
                        print(ind3_md)
                        inds3_final = (ind3_mn + ind3_md) / 2
                        print(inds3_final)
                        inds3.append(inds3_final) 

                    else:
                        if cls[2] < correct_cluster_ground_truth[2]:
                            ind3_mn = round(cls[2] / correct_cluster_ground_truth[2], 3) 
                            if ind3_mn > 1:
                                ind3_mn = 1
                        else:
                            ind3_mn = round(correct_cluster_ground_truth[2] / cls[2], 3) 
                            if ind3_mn > 1:
                                ind3_mn = 1
                        print(ind3_mn)
                        if cls[3] < correct_cluster_ground_truth[3]:
                            ind3_md = round(cls[3]/correct_cluster_ground_truth[3], 3) 
                            if ind3_md > 1:
                                ind3_md = 1
                        else:
                            ind3_md = round(correct_cluster_ground_truth[3]/cls[3], 3) 
                            if ind3_md > 1:
                                ind3_md = 1
                        print(ind3_md)
                        inds3_final = (ind3_mn + ind3_md) / 2
                        print(inds3_final)
                        inds3.append(inds3_final) 

                
                indicator2 = round((sum(inds2) / len(inds2)), 3)
                print("final indicator2 for polygon: ")
                print(indicator2)
                
                indicator3 = round((sum(inds3) / len(inds3)),3)
                print("final indicator3 for polygon: ")
                print(indicator3)
        polygons_indicators.append([p[0], p[-1], indicator2, indicator3])
    return polygons_indicators


#clusters_num(pls_case1, ground_truth_case1, clust_approach2_case1)
lists_lines = indicators(pls_case, ground_truth_case2, clust_approach)
polygon_gtlines = lists_lines[0]
polygon_clustlines = lists_lines[1]

print("Start comparing for indicator 2 And indicator 3!!")
a = indicators_comparison(polygon_gtlines, polygon_clustlines)

in1 = 0
in2 = 0
in3= 0
ind1_l = []
ind1_h = []
ind2_l = []
ind2_h = []
ind3_l = []
ind3_h = []

for aaa in a:
    print(aaa)

    in1 += aaa[1]
    if aaa[1] < 0.5:
        ind1_l.append(1)
    elif aaa[1] > 0.9:
        ind1_h.append(1)

    in2 += aaa[2]
    if aaa[2] < 0.5:
        ind2_l.append(1)
    elif aaa[2] > 0.9:
        ind2_h.append(1)

    in3 += aaa[3]
    if aaa[3] < 0.5:
        ind3_l.append(1)
    elif aaa[3] > 0.9:
        ind3_h.append(1)

print(round(in1/5, 2))
print(round(in2/5, 2))
print(round(in3/5, 2))

final_idn = round((round(in1/5, 3) * 0.25 + round(in2/5, 3) * 0.5 + round(in3/5, 3) * 0.25), 2)
print(final_idn)

print(len(ind1_l))
print(len(ind1_h))
print(len(ind2_l))
print(len(ind2_h))
print(len(ind3_l))
print(len(ind3_h))

print("check1")

tin = 0
for aa in a:
    print(aa)
    print(aa[0])
    index = aa[1]*0.25 + aa[2]*0.5 + aa[3]*0.25
    print(index)
 
    tin+= index

print("final index for polygons of that case: ")
print(tin)
print(tin/5)

# check relation between length and number of clusters
def checkrelation(pls,edges,clustlines):
    l_1 = []
    l_2 = []
    l_3 = []
    l_4 = []
    l_o4 = []
    for p in range(len(pls)):
        rec = pls.loc[p, 'geometry']
        print(rec)
        n_ed = edges.sindex.query(rec, predicate = 'intersects')
        n_cl = clustlines.sindex.query(rec, predicate = 'intersects')
        if len(n_ed) != 1:
            print("EEEE MLKIA EEEEEE !!!! OOPP!!!!!!")
        #print("this are the edges: ")
        for i in n_ed:
            r = edges.loc[i,'geometry']
            #print(r)
        #print("this are the clust lines top: ")
        for ii in n_cl:
            rr = clustlines.loc[ii, 'geometry']
            #print(rr)
        

        print(len(n_cl))
        if len(n_cl) == 1:
            l_1.append(edges.loc[n_ed[0],'geometry'].length)
        elif len(n_cl) == 2:
            l_2.append(edges.loc[n_ed[0],'geometry'].length)
        elif len(n_cl) >2 and len(n_cl) <=5:
            l_3.append(edges.loc[n_ed[0],'geometry'].length)
        elif len(n_cl) > 5 and len(n_cl) <=10:
            l_4.append(edges.loc[n_ed[0],'geometry'].length)
        elif len(n_cl) > 10:
            l_o4.append(edges.loc[n_ed[0],'geometry'].length)

    a = areameanmedian(l_1)
    m_1 = a[0]
    md_1 = a[1]
    print(m_1)
    print(md_1)
    aa = areameanmedian(l_2)
    m_2 = aa[0]
    md_2 = aa[1]
    print(m_2)
    print(md_2) 
    aaa =  areameanmedian(l_3)
    m_3 = aaa[0]
    md_3 = aaa[1]
    print(m_3)
    print(md_3)
    aaaa = areameanmedian(l_4)
    m_4 = aaaa[0]
    md_4 = aaaa[1]
    print(m_4)
    print(md_4)
    aaaaa = areameanmedian(l_o4)
    m_5 = aaaaa[0]
    md_5 = aaaaa[1]
    print(m_5)
    print(md_5)


# not used
def stats_ofnewlines(pls, ground_truth, clust_approach):
    polygon_gtlines = []
    polygon_clustlines = []
    for i in range(len(pls)):
        truth = []
        clust = []
        rec = pls.loc[i, 'geometry']
        #print("For polygon: ", rec)
        n = ground_truth.sindex.query(rec, 'intersects')
        n2 = clust_approach.sindex.query(rec, 'intersects')
        for ind in n:
            ln = ground_truth.loc[ind, 'geometry']
            truth.append(ln)
        for ind2 in n2:
            ln2 = clust_approach.loc[ind2, 'geometry']
            clust.append(ln2)
        #Find the lenght of the lines for ground truth and for clustering approach for each polygon
        lengths_of_lines_t = []
        idd = 0
        for ll in truth:
            #print(ll)
            lengg = round(ll.length, 2)
            how_many_measure = round(lengg)
            interval = 1
            measure_lines_width = []
            for hh in range(how_many_measure):
                m = LineString([ll.interpolate(interval-1, normalized = False), ll.interpolate(interval, normalized = False)])
                #print(m)
                interval+=1
                mdpoint = m.interpolate(0.5, normalized = True)
                leftline = m.parallel_offset(50, 'left')
                rightline = m.parallel_offset(50, 'right')
                midpointleft = leftline.interpolate(0.5, normalized = True)
                midpointright = rightline.interpolate(0.5, normalized = True)
                prependicular_line = LineString([midpointleft, midpointright])
                prep_line = [prependicular_line, 1]
                #print(prep_line[0])
                #for each prep line find the 2 closest points to the polygon that contains the initial line (create the measuring line)
                inn = pls.sindex.query(mdpoint, predicate="within")
                if len(inn)!= 1 and len(inn)!= 0:
                    print("error")
                if len(inn)==0:
                    continue
                pol_contain_line = pls.loc[inn[0], 'geometry']
                oww = width(prep_line[0], pol_contain_line, ll, prep_line[1])
                measure_lines_width.append(oww[0])
            s = areameanmedian(measure_lines_width)
            lengths_of_lines_t.append([lengg, round(s[0],3), round(s[1],3), idd])
            idd+=1
        #print("Ground truth results to: ", len(truth), " lines, with lengths, mean width and median width: ")
        polygon_gtlines.append([rec, len(truth), lengths_of_lines_t])
        for op in lengths_of_lines_t:
            print(op)
            

        lengths_of_lines_c = []
        idd2 = 0
        for ll2 in clust:
            #print(ll)
            lengg2 = round(ll2.length, 2)
            how_many_measure2 = round(lengg2)
            interval2 = 1
            measure_lines_width2 = []
            for hh2 in range(how_many_measure2):
                m2 = LineString([ll2.interpolate(interval2-1, normalized = False), ll2.interpolate(interval2, normalized = False)])
                #print(m)
                interval2+=1
                mdpoint2 = m2.interpolate(0.5, normalized = True)
                leftline2 = m2.parallel_offset(50, 'left')
                rightline2 = m2.parallel_offset(50, 'right')
                midpointleft2 = leftline2.interpolate(0.5, normalized = True)
                midpointright2 = rightline2.interpolate(0.5, normalized = True)
                prependicular_line2 = LineString([midpointleft2, midpointright2])
                prep_line2 = [prependicular_line2, 1]
                #print(prep_line[0])
                #for each prep line find the 2 closest points to the polygon that contains the initial line (create the measuring line)
                inn2 = pls.sindex.query(mdpoint2, predicate="within")
                if len(inn2)!= 1 and len(inn2)!= 0:
                    print("error")
                if len(inn2)==0:
                    continue

                pol_contain_line2 = pls.loc[inn2[0], 'geometry']
                oww2 = width(prep_line2[0], pol_contain_line2, ll2, prep_line2[1])
                measure_lines_width2.append(oww2[0])
            s2 = areameanmedian(measure_lines_width2)
            lengths_of_lines_c.append([lengg2, round(s2[0],3), round(s2[1],3), idd2])
            idd2+=1
        polygon_clustlines.append([rec, len(clust), lengths_of_lines_c])
        #print("Approach results to: ", len(clust), " lines, with lengths, mean width and median width: ")
        for op2 in lengths_of_lines_c:
            #print(op2)
            continue    

    return polygon_gtlines, polygon_clustlines   
 
def findsimilarities(polygon_gtlines, polygon_clustlines):
    identical = []
    quite_similar = []
    no_similar = []
    for p in polygon_gtlines:
        print(p)
        for pp in polygon_clustlines:
            #find the same polygon for gt and for clust approach
            if p[0] == pp[0]:
                print(pp)
                print(p[0])
                similarity1 = []
                similarity2 = []
                print("num of clusters for gt: ", p[1])
                print("num of clusters for clustering approach: ", pp[1])
                # if we have the same num of clusters
                if p[1] == pp[1]:
                    print(p[1])
                    for ii in p[2]:
                        for iii in pp[2]:
                            #check that we are comparing the correct lines based on the id
                            if ii[3] == iii[3]:
                                print(ii)
                                print(iii)
                            #if length of lines is almost similar less than 5 m and mean/median width is almost similar 
                                if abs(ii[0]-iii[0]) < 5:
                                    if abs(ii[1]-iii[1]) < 0.5 and abs(ii[1]-iii[1]) < 0.5:
                                        print("the line is identical")
                                        similarity1.append(1)
                                else:
                                    continue
                    if len(similarity1) == p[1]:
                        print(p[1], "lines are identical so identical")
                        identical.append(p[0])
                    elif len(similarity1) == (p[1] -1) :
                        print(len(similarity1), "lines are identical so quite similar")
                        quite_similar.append(p[0])
                    else:
                        print(len(similarity1), "lines are identical so unsimilar")
                        no_similar.append(p[0])
                
                #if we have +- num of clusters
                elif abs(p[1]-pp[1]) == 1:
                    print("no same num of clusters")
                    for ck in p[2]:
                        for ck2 in pp[2]:
                            print(ck[0])
                            print(ck2[0])
                            if abs((ck[0] - ck2[0])) < 6 and abs((ck[1]- ck2[1])) < 1 and abs((ck[2]-ck2[2])<1):
                                print("they are quite similar")
                                if ck not in similarity2:
                                    similarity2.append(ck)

                    if len(similarity2) == p[1]:
                        print(len(similarity2), "lines are quite similar so they are similar")
                        quite_similar.append(p[0])
                    else:
                        print(len(similarity2), "lines are quite similar so they are unsimilar")
                        no_similar.append(p[0])
                else:
                    print("they are unsimilar")
                    no_similar.append(p[0])
    return [identical, quite_similar, no_similar]



