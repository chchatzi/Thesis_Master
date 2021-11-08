from re import M
import csv
import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#from matplotlib import figure
from shapely.geometry import mapping, Polygon
import fiona
from geomet import wkt
import pandas as pd
import geopandas as gpd
from geopandas import GeoDataFrame
import shapefile
from shapely.geometry import box
from shapely.validation import explain_validity
import networkx as nx
import osmnx as ox
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.geometry import MultiLineString
from shapely.geometry import Point
from operator import itemgetter
from shapely.ops import linemerge
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
import os
from shapely.ops import linemerge
import statistics
from shapefile_manipulation import df_roads, epgs, df_cross, df_another, df_rounds, df_t, df_x, dist_for_offsets

'''
Width clustering and width estimations performed!

Input --> Takes the final road edge polygons as they result from shapefile_manipulation.py file
Output --> Clustered centerlines based on width measurements, assigned with width statistics (mean,median, max,min etc.)

Predefined values:

-interval =  Corresponds to the measuring interval that will be used for width estimation. (every how m the centerlines will be 'cut')
'''

#Help functions
def removedoubles(pols_shapefile):
    all_p = []
    for i in range(len(pols_shapefile)):
        rec = pols_shapefile.loc[i, 'geometry']
        all_p.append(rec)

    for i in all_p:
        if all_p.count(i) > 1:
            all_p.remove(i)

    final_p = list_to_dgf(all_p,epgs)
    final_p.to_file("final_polygons_for_clstering.shp")

def fivelines(line):
    tenlines = []
    line1 = [LineString([line.interpolate(0, normalized = True), line.interpolate(0.2, normalized = True)]), 1]
    line2 = [LineString([line.interpolate(0.2, normalized = True), line.interpolate(0.4, normalized = True)]), 2]
    line3 = [LineString([line.interpolate(0.4, normalized = True), line.interpolate(0.6, normalized = True)]), 3]
    line4 = [LineString([line.interpolate(0.6, normalized = True), line.interpolate(0.8, normalized = True)]), 4]
    line5 = [LineString([line.interpolate(0.8, normalized = True), line.interpolate(1, normalized = True)]), 5]
    tenlines.append(line1)
    tenlines.append(line2)
    tenlines.append(line3)
    tenlines.append(line4)
    tenlines.append(line5)
    return tenlines
    
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
    paro = len(list_areas)-1
    std = round(math.sqrt(s/paro), 4)
    #print("this is std: ", std)
    return [mean, median, std]

def generatestats2(l,nm):
    list_lengths = []
    smalls = []
    bigs = []
    for i in range(len(l)):
        ln = l.loc[i, 'geometry']
        list_lengths.append(ln.length)
    
    a = areameanmedian(list_lengths)
    mean = a[0]
    median = a[1]
    std = a[2]
    print(mean -std)
    print(mean+std)
    for i in range(len(l)):
        ln = l.loc[i, 'geometry']
        if ln.length < mean - (mean*0.9):
            smalls.append(ln)
        if ln.length > 2*mean:
            bigs.append(ln)
        
    print(len(smalls))
    print(len(bigs))
    p_small = (len(smalls)*100) / len(l) 
    p_big = (len(bigs)*100) / len(l)
    max_l = max(list_lengths)
    min_l = min(list_lengths)

    with open(nm, 'w', newline='') as fileeee:
        fieldnames = ['total_lines', 'Mean length (m)', 'Std', 'Median length (m)','Max centerline',
         'Min centerline', 'Total small centerlines', '%small centerlines', 'Total big lines', '%big lines']
        writer = csv.DictWriter(fileeee, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'total_lines':len(l), 'Mean length (m)': mean, 'Std':std, 'Median length (m)':median,'Max centerline':max_l,
        'Min centerline':min_l, 'Total small centerlines':len(smalls), '%small centerlines':p_small, 'Total big lines':len(bigs), '%big lines':p_big})

def avgwidth2(lst):
    s = 0
    for i in lst:
        s += i
    avg = s / len(lst)
    return round(avg, 4)

def med_min_max_range_mean_std(lst):
    s = 0
    for i in lst:
        s+=i
    mean = s/len(lst)

    s = 0
    for l in lst:
        d = (l - mean) ** 2
        s += d
    paro = len(lst)-1
    std = round(math.sqrt(s/paro), 4)
    med = statistics.median(lst)
    maxi = max(lst)
    mini = min(lst)
    rangee = maxi-mini
    #print("this is median: ", med)
    #print("this is min: ", mini)
    #print("this is max: ", maxi)
    return med,mini,maxi,rangee,mean,std

def list_to_dgf(listt,epgs):
    d = {'geometry':[]}
    for pp in listt:
        d['geometry'].append(pp)

    gdf = gpd.GeoDataFrame(d, crs = epgs)
    return gdf

def delete_files_111():
    os.remove("only_roads_after_recreation.shp")
    os.remove("only_roads_after_recreation.cpg")
    os.remove("only_roads_after_recreation.dbf")
    os.remove("only_roads_after_recreation.shx")

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

def stats_of_measure_lines(all_widths):
    t = 0
    for w in all_widths:
        t+=w
    
    mm = t / len(all_widths)
    s = 0
    for ww in all_widths:
        d = (ww - mm) ** 2
        s += d
    paro = len(all_widths)-1
    std = round(math.sqrt(s/paro), 4)
    print("This is the mean lenght of all measuring lines : ", mm)
    print("This is the std from the mean length of measure lines: ", std)
    return [mm, std]


#Main functions
def findmeasurelines(lns, pls, dist,epgs, interval):
    all_widths_of_dataset = []
    all_m = [] #contains all the measuring lines for all the centerlines of the dataset (width measurement, wkt)
    points = []# midpoints of each measuring lines for all centerlines of dataset ([pointx,pointy, width measurement])
    all_small_centerlines = [] #contains all the parts of the centerlines after they cut based on a measuring interval
    doubles = []
    print("we have", len(lns), "centerlines and, ", len(pls), "road edge polygons")
    #for each centerline in lns geodataframe, find the number of measuring lines (every x meters) based on their length
    for i in range(len(lns)):
        rec = lns.loc[i, 'geometry']
        #print("this is centerline: ")
        #print(rec)
        lengg = rec.length
        how_many_measure = round(lengg/interval)
        #print("lenght of centerline: ", lengg, " number of measuring lines for this centerline: ", how_many_measure)       
        nm = interval
        id = 0
        midpoints = []
        #print("measuring lines: ")
        if how_many_measure < 1:
            print("super small centerline!")
            print(rec)
            continue
        else:
            for ii in range(how_many_measure):
                #cut each line every x(interval) m to smaller parts
                m = [LineString([rec.interpolate(nm-interval, normalized = False), rec.interpolate(nm, normalized = False)]), id]
                nm+=interval
                id+=1
                mdpoint = m[0].interpolate(0.5, normalized = True)
                #everty x m get a prependicular line
                leftline = m[0].parallel_offset(dist, 'left')
                rightline = m[0].parallel_offset(dist, 'right')
                midpointleft = leftline.interpolate(0.5, normalized = True)
                midpointright = rightline.interpolate(0.5, normalized = True)
                prependicular_line = LineString([midpointleft, midpointright])
                prep_line = [prependicular_line, m[1]]
                #print(prep_line[0])
                #for each prep line find the 2 closest points to the polygon that contains the initial line (create the measuring line)
                inn = pls.sindex.query(mdpoint, predicate="within")
                #print(len(inn))
                if len(inn) < 1:
                    #in some cases the first/last measuring line is not inside any polygon (centerline is not starting/ending exactly where the polygon starts ends)
                    #thus in that cases we continue with the next measuring line
                    continue
                elif len(inn) > 1:
                    #in cases of flyovers, some polygons are overlapping. Our methodology is not developed to work for those cases but only for the ordinary roads
                    continue
                else:
                    #find the polygon that contains the part of centerline (as it has been created previously)
                    pol_contain_line = pls.loc[inn[0], 'geometry']
                    #find the measuring line (intersection of prependicular line with polygon that contains the line)
                    oww = width(prep_line[0], pol_contain_line, rec, prep_line[1])
                    #check if oow[1] (the measuring line) intersects with more than 1 centerlines
                    nnn = lns.sindex.query(oww[1],  predicate="intersects")
                    if len(nnn)>1:
                        #print(oww[1])
                        doubles.append([oww[1], m[0]])
                    if oww[1].geom_type != 'LineString':
                        for iii in oww[1]:
                            if iii.geom_type == 'LineString':
                                midpoint = iii.interpolate(0.5, normalized = True)
                                midpoints.append([midpoint, m[1]])
                                point = [midpoint.x, midpoint.y, oww[0]]
                                all_m.append(oww)
                                points.append(point)
                                all_widths_of_dataset.append(oww[0])
                                all_small_centerlines.append(m[0])
                    else:
                        midpoint = oww[1].interpolate(0.5, normalized = True)
                        #print(oww[1])
                        #print(midpoint)
                        #print(m[1])
                        midpoints.append([midpoint, m[1]])
                        point = [midpoint.x, midpoint.y, oww[0]]
                        all_m.append(oww)
                        points.append(point)
                        all_widths_of_dataset.append(oww[0])
                        all_small_centerlines.append(m[0])
    
        #check the average euclidean distance (x,y) between midpoints of succesive measuring lines for each centerline--> use it for distance threshold for clustering
        for pntt in range(len(midpoints)-2):
            pp1 = midpoints[pntt]
            pp2 = midpoints[pntt+1]
            if abs(pp1[1] - pp2[1]) == 1:
                p1 = np.array((pp1[0].x, pp1[0].y))
                p2 = np.array((pp2[0].x, pp2[0].y))
                ddd = np.linalg.norm(p1-p2)
                if ddd < interval + 0.1 and ddd > interval - 0.1 :
                    point1 = [pp1[0].x, pp1[0].y]
                    point2 = [pp2[0].x, pp2[0].y]
                    typ_d = np.linalg.norm(np.array(point1)-np.array(point2))

    
    #based on all the measuring lines of the dataset compute some statistics ---> to generate the distance threshold for clustering
    aou = stats_of_measure_lines(all_widths_of_dataset)
    mdmnmx = med_min_max_range_mean_std(all_widths_of_dataset)
    mm = aou[0]
    mdd = mdmnmx[0]
    mnn = mdmnmx[1]
    mxx = mdmnmx[2]
    std = aou[1]
    help_lst = []
    for jk in all_m:
        if jk[1].geom_type != 'LineString':
            print("i have a measuring line that is a Geometry collection. Propably belongs to an intersection polygon that is not identified correctly!")
            for iiii in jk[1]:
                if iiii.geom_type == 'LineString':
                    help_lst.append(iiii)
        else:
            help_lst.append(jk[1])
    all_m_g = list_to_dgf(help_lst, epgs)
    small_centerlines = list_to_dgf(all_small_centerlines,epgs)

    return points, all_m_g, point1, point2, std, typ_d, small_centerlines, doubles


def find_distance_threshold_fixedwidth(p1, p2, width_value, typ_d):
    '''Based on the typical distance of 2 succesive measuring lines (results based on the measuring interval + a small tolerance)
    --> the distance threshold for clustering is generated.'''
    p1_n = np.array((p1[0], p1[1], 0))
    p2_n = np.array((p2[0], p2[1], width_value))
    ddd = np.linalg.norm(p1_n-p2_n)
    f_d = ddd-0.3
    print("this is distance between 2 typical points in euclidean space : ", typ_d )
    print("this is the width value that will be used as width threshold for clusters: ", width_value)
    print("this is the distance threshold for clustering: ", f_d)
    return f_d


def clusterpoints(points, lines_gdf, epgs, thres, s_ln, doubles):
    ar = np.array(points)
    new_lines = []
    print("Start clustering")
    print("Midpoints of measuring lines that are going to be used for clustering (x,y,width): ")
    print(ar)
    print("number of points: ", len(ar))
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(ar[:,0], ar[:,1], ar[:,2], s=50)
    ax.view_init(azim=200)

    linked = linkage(ar, 'single')
    plt.figure(figsize=(10, 7))
    dendrogram(linked,
                orientation='top',
                distance_sort='descending',
                show_leaf_counts=True)
    #plt.show()
    cluster = AgglomerativeClustering(n_clusters = None, distance_threshold= thres, affinity='euclidean', linkage='single')
    y_hc = cluster.fit_predict(ar)
    #print(y_hc)
    cluster_index = []
    for pnt in y_hc:
        if pnt not in cluster_index:
            cluster_index.append(pnt)
    print("this is the number of clusters: ", len(cluster_index))
    lines_gdf.to_file("clustering_measuring_lines.shp")

    #add clustered points of the same index to the same cluster (list = points_of_cluster_x)
    f = dict()
    for i in cluster_index:
        f[i] = []
        #print("NEO CLUSTER : ")
        points_of_cluster_x = []
        points_of_clust = ar[y_hc == i]
        for p in points_of_clust:
            pnt = Point(p[0],p[1])
            #print(pnt)
            points_of_cluster_x.append(pnt)
        pnt_gdf = list_to_dgf(points_of_cluster_x,epgs)
        #for each clustered mid point find the corresponding measuring line (from lines_gdf)
        for ox in range(len(pnt_gdf)):
            recc = pnt_gdf.loc[ox, 'geometry']
            l_ind = lines_gdf.sindex.query(recc)
            #print(len(l_ind))
            #if the midpoint intersects with more than one measuring lines , find the correct one by using the distance of midpoint with the midpoint of each measure line (with the correct one the distance should be 0)
            if len(l_ind) == 2:
                #print(len(l_ind))
                a1 = lines_gdf.loc[l_ind[0], 'geometry']
                a2 = lines_gdf.loc[l_ind[1], 'geometry']
                #print(a1)
                #print(a2)
                #print(recc)                
                md_p1 = a1.interpolate(0.5, normalized = True)
                md_p2 = a2.interpolate(0.5, normalized = True)
                d1 = recc.distance(md_p1)
                d2 =  recc.distance(md_p2)
                #print(d1)
                #print(d2)
                if d1 < d2:
                    f[i].append(lines_gdf.loc[l_ind[0], 'geometry'])
                    #print(lines_gdf.loc[l_ind[0], 'geometry'])
                else:
                    f[i].append(lines_gdf.loc[l_ind[1], 'geometry'])
                    #print(lines_gdf.loc[l_ind[1], 'geometry'])

            elif len(l_ind) == 3:
                #print(len(l_ind))
                a1 = lines_gdf.loc[l_ind[0], 'geometry']
                a2 = lines_gdf.loc[l_ind[1], 'geometry']
                a3 = lines_gdf.loc[l_ind[2], 'geometry']
                #print(a1)
                #print(a2)
                #print(a3)
                #print(recc)                
                md_p1 = a1.interpolate(0.5, normalized = True)
                md_p2 = a2.interpolate(0.5, normalized = True)
                md_p3 = a3.interpolate(0.5, normalized = True)
                d1 = recc.distance(md_p1)
                d2 =  recc.distance(md_p2)
                d3 = recc.distance(md_p3)
                #print(d1)
                #print(d2)
                #print(d3)
                if d1 < d2 and d1 < d3:
                    f[i].append(lines_gdf.loc[l_ind[0], 'geometry'])
                    #print(lines_gdf.loc[l_ind[0], 'geometry'])
                elif d2<d1 and d2<d3:
                    f[i].append(lines_gdf.loc[l_ind[1], 'geometry'])
                    #print(lines_gdf.loc[l_ind[1], 'geometry'])
                else:
                    f[i].append(lines_gdf.loc[l_ind[2], 'geometry'])
                    #print(lines_gdf.loc[l_ind[2], 'geometry'])
            else:
                #print("ok")
                #print(lines_gdf.loc[l_ind[0], 'geometry'])
                f[i].append(lines_gdf.loc[l_ind[0], 'geometry'])
    print(f)
    print("GO")
    #from clustered measuring lines find the corresponding small centerlines.
    #if a measuring line intersects with more than one small centerline (use doubles list to find the correct one)
    final = []
    idd = -1
    for clstr in f:
        idd+= 1
        print("clusterr")
        cluster_lines = []
        cluster_lengths = []
        for ln in f[clstr]:
            print("measuring line: ")
            print(ln)
            n = s_ln.sindex.query(ln,  predicate="intersects")
            if len(n) != 1:
                #print("clusteredd measuring line intersects with more than one small centeline part: ")
                #print(len(n))
                #print(ln)
                for db in doubles:
                    if ln == db[0]:
                        small_line = db[1]
                        print("correspond to small centerline: ")
                        print(small_line)
                        cluster_lines.append(small_line)
                        cluster_lengths.append(ln.length)
                    else:
                        continue
            else:
                small_line = s_ln.loc[n[0], 'geometry']
                print("correspond to small centerline: ")
                print(small_line)
                cluster_lines.append(small_line)
                cluster_lengths.append(ln.length)
        
        if len(cluster_lines) == 0:
            print("i am here")
            continue
        if len(cluster_lines) == 1:
            #print("only one line in the cluster")
            #print(cluster_lengths[0])
            m = cluster_lines[0]
            men = cluster_lengths[0]
            med = "only one line in the cluster"
            mn = "only one line in the cluster"
            mx = "only one line in the cluster"
            rang = "only one line in the cluster"
            stdd = 'only one line in the cluster'

        else:
            m_mn_mx =  med_min_max_range_mean_std(cluster_lengths)
            men = m_mn_mx[4]
            rang = m_mn_mx[3]
            stdd = m_mn_mx[5]
            med = m_mn_mx[0]
            mn = m_mn_mx[1]
            mx = m_mn_mx[2]
            #print("final")
            m = linemerge(cluster_lines)
            print(m)

        final.append([m, men, med, mn, mx, idd, rang,stdd])

    return final


def writeclusterlines2(final_list):
    alllins = final_list
    schema2 = {
        'geometry': 'LineString',
        'properties': {'id' : 'int', 'Mean':'float:13.3','STD':'str', 'Median':'str','Min':'str', 'Max':'str','Range':'str'},
    }

    with fiona.open('clustered_centerlines.shp', 'w', 'ESRI Shapefile', schema2) as c:
        ## If there are multiple geometries, put the "for" loop here
        for ii in alllins:
            c.write({
                'geometry': mapping(ii[0]),
                'properties': {'id' : ii[5], 'Mean':ii[1], 'STD':str(ii[7]) , 'Median':str(ii[2]),'Min':str(ii[3]), 'Max':str(ii[4]), 'Range': str(ii[6])},
            })


interval = 5
pols_shapefile = df_roads
lns_shapefile = "edges_projected.shp"

#Just in case there was an error during previous process of re-create road polygons --> check if there is dublicates an remove
#Not nessecary step maybe need to be removed (it slows down a lot the process for big dataset) + No case of dublicates for the tested areas 
#print("Remove dublicate polygons")
#removedoubles(pols_shapefile)
#
print("Each centerline will be 'cut' (for width measurements) based on a measuring interval of : ", interval, "meters!")

pols_shapefile3 = "final_polygons_for_clstering.shp" 

pls = gpd.read_file(pols_shapefile3)
lns = gpd.read_file(lns_shapefile)
print("START CLUSTERING!!")
#find the measuring lines for all the dataset
print("Start find the measuring lines for all the centerlines of the dataset at equal intervals")
pl = findmeasurelines(lns,pls, dist_for_offsets, epgs, interval)
point1 = pl[2]
point2 = pl[3]
typ_d = pl[5]
#Generate the distance threshold for clustering --> based on measuring interval + predefined width value (2.0 used)
thres = find_distance_threshold_fixedwidth(point1, point2, 2.0, typ_d)
#perform the clustering
c = clusterpoints(pl[0],pl[1], epgs, thres , pl[6], pl[7])


#write the file
writeclusterlines2(c)


#delete_files_111()
e = gpd.read_file('clustered_centerlines.shp')
generatestats2(e,'clustering_lines_stats.csv')

#not used
def find_distance_threshold(p1, p2, std, typ_d):    

    p1_n = np.array((p1[0], p1[1], 1))
    p2_n = np.array((p2[0], p2[1], 1+std/2))
    ddd = np.linalg.norm(p1_n-p2_n)
    print("this is distance between 2 typical points in euclidean space : ", typ_d )
    print("this is the std/2 that will be used as width threshold for clusters: ", std/2)
    print("this is the distance threshold for clustering: ", ddd)
    return ddd+0.2

def findextremepoints(points):
    #for ii in points:
        #print(ii)
    if len(points) == 1:
        return 1
    elif len(points) == 2:
        #ln = LineString([points[0], points[1]])
        #print(ln)
        return 2
    else:
        distancess = []
        p1 = points[0]
        for i in range(len(points)):
            d = p1.distance(points[i])
            distancess.append([d, points[i]])
        sorted_dist = sorted(distancess, key=lambda x: x[0])
        extreme_p = sorted_dist[-1][1]
        #print("extreme poitns")
        #print(extreme_p)
        distances2 = []
        for io in range(len(points)):
            d = extreme_p.distance(points[io])
            distances2.append([d, points[io]])
        sorted_dist2 = sorted(distances2, key=lambda x: x[0])
        extreme_p2 = sorted_dist2[-1][1]
        linee = [Point(extreme_p.x, extreme_p.y)]
        for pp in sorted_dist2:
            p = Point(pp[1].x, pp[1].y)
            linee.append(p)
        #print(extreme_p2)
        #print(linee)
        ln = LineString(linee)
        #print(ln)
        return ln

def merge_lines_merge_polygons(pols_file, epgs):
    # merge polygons
    p = gpd.read_file(pols_file)
    p['d_field'] = 1
    p_new = p.dissolve(by='d_field')
    p_new.crs = epgs
    p_new.to_file("diss_polygons.shp")
    nn = p_new.to_crs("EPSG:4326")
    nn.to_file("wgs84_dissolved_polygons.shp")
    
    d_polygons = gpd.read_file("wgs84_dissolved_polygons.shp")
    pl = d_polygons.loc[0,'geometry']
    #print(pl.bounds)
    b = box(pl.bounds[0],pl.bounds[1],pl.bounds[2],pl.bounds[3])
    #print(b)
    #Generate the graph WGS84
    G = ox.graph_from_polygon(b, network_type='drive')
    #ox.plot_graph(G)
    ox.io.save_graph_shapefile(G,"graph")
    # merge lines
    aa = gpd.read_file("C:/Users/chcha/Desktop/GEOMATICS/2nd Year/Thesis/p2_technical/code/methodology1/graph/edges.shp")
    lines = aa.to_crs(epgs)
    #print(lines.crs)
    lstt = []
    for i in range(len(lines)):
        rec = lines.loc[i]
        if rec['highway'] == 'motorway_link' or rec['highway'] == 'motorway' or rec['highway'] == 'trunk_link' or rec['highway'] == 'trunk':
            continue
        else:
            nm = rec['name']
            g = rec['geometry']
            lstt.append([nm,g])
    f = dict()
    for opa in lstt:
        #print(opa)
        if opa[0] not in f:
            f[opa[0]] = [opa[1]]
        else:
            f[opa[0]].append(opa[1])
    merged_lines = []
    for ii in f:
        m = linemerge(f[ii])
        merged_lines.append(m)
    vcc = list_to_dgf(merged_lines, epgs)
    vcc.to_file("clusterring_merged_lines.shp")

    return "clusterring_merged_lines.shp", vcc, "diss_polygons.shp", p_new