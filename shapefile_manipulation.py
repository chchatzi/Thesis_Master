from numpy import min_scalar_type
import shapefile
import shapely.wkt
import math
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.geometry import MultiLineString
from shapely.geometry import Point
from shapely.validation import explain_validity
import matplotlib.pyplot as plt
import geopandas as gpd
import os
import statistics
import csv
from recreate_polygons import recreate_int_polygons, list_to_dgf, close_holes
from merge_lines_based_nameandhighway import merge_lines_namehighway

'''
INPUT --> Polygons shapefile + some predifined thresholds (see below)
Returns --> Creates 6 new polygon shapefiles. 5 shapefiles with the different intersection types (X,Cross,Roundabouts,T,Another_type) + 1 shapefile with final road edge polygons (motorways removed)
Returns --> edges projected, original road centerlines (downloaded from OSM) for the selected area (final edges that will be used for clustering and width estimation)


Recreate road polygons --> Uses recreate_int_polygons function from recreate_polygons.py file
Identify intersections --> X,Cross,Roundabouts,T, every other intersection --> another_type

Predefined values needed for polygons re-creation:

- epgs = EPSG coordinate reference system (i.e for Helsinki = "EPSG:3879")

-bufer =  Buffer size for polygon re-creation (19-21 used)

-dist_for_split = Distance threshold for spliting the large re-created intersection polygons that contains 2 nodes into 2 smaller int polygons (25-30 ideal). 

-area_for_mrg = Area threshold for define super small polygons that result from difference between polygons and buffer. If a polygon small that thres then merge it with
a neighboring polygon. (30-40 used)

-dist_for_offsets = The distance that the offest for defining the meauring lines. (Any number over 10m will be fine. Just need to make sure that offsets will lie outside of
the polygon)

-polygons_file =  Shapefile that contains the polygons that is to be re-created, and identify intersections
'''

#Intersection identification functions

def checkforTcase2(twentylines, roadline, size):
    if size == 2: 
        if len(twentylines) > 19:
            first = [twentylines[0][0], twentylines[1][0], twentylines[2][0]]#, twentylines[3], twentylines[4], twentylines[5]]
            mid = [twentylines[8][0], twentylines[9][0], twentylines[10][0], twentylines[11][0]]
            last = [twentylines[17][0], twentylines[18][0], twentylines[19][0]]
        elif len(twentylines) > 18:
            first = [twentylines[0][0], twentylines[1][0], twentylines[2][0]] #, twentylines[3], twentylines[4], twentylines[5]]
            mid = [twentylines[8][0], twentylines[9][0], twentylines[10][0], twentylines[11][0]]
            last = [twentylines[16][0], twentylines[17][0], twentylines[18][0]]
        elif len(twentylines) > 17:
            first = [twentylines[0][0], twentylines[1][0], twentylines[2][0]] #, twentylines[3], twentylines[4], twentylines[5]]
            mid = [twentylines[7][0], twentylines[8][0], twentylines[9][0], twentylines[10][0], twentylines[11][0]]
            last = [twentylines[15][0], twentylines[16][0], twentylines[17][0]]
        else:
            return False
    elif size == 1:
        if len(twentylines) > 9:
            first = [twentylines[0][0], twentylines[1][0], twentylines[2][0]]
            mid = [twentylines[5][0], twentylines[6][0]]
            last = [twentylines[8][0], twentylines[9][0]]
        elif len(twentylines)> 8:
            first = [twentylines[0][0], twentylines[1][0]]
            mid = [twentylines[4][0], twentylines[5][0]]
            last = [twentylines[7][0], twentylines[8][0]]
        elif len(twentylines)> 7:
            first = [twentylines[0][0], twentylines[1][0]]
            mid = [twentylines[4][0], twentylines[5][0]]
            last = [twentylines[6][0], twentylines[7][0]]
        else:
            return False

    f = avgwidth2(first)
    m = avgwidth2(mid)
    l = avgwidth2(last)
    #print(f, "",m, "", l)
    if m>f*1.5 and m>l*1.5:
        #print("mid is 1.5 times higher")
        #check if it is a T intersection
        #take 4 measuring lines that are in the midlle
        #take the 2 end points of each line
        if size==2:
            end8one = Point(twentylines[8][1].coords[0])
            end8tow = Point(twentylines[8][1].coords[-1])
            end9one = Point(twentylines[9][1].coords[0])
            end9tow = Point(twentylines[9][1].coords[-1])
            end10one = Point(twentylines[10][1].coords[0])
            end10tow = Point(twentylines[10][1].coords[-1])
            end11one = Point(twentylines[11][1].coords[0])
            end11tow = Point(twentylines[11][1].coords[-1])

            #find the midpoints of the road lines that are in line with the 4 mid lines
            line9 = LineString([roadline.interpolate(0.4, normalized = True), roadline.interpolate(0.45, normalized = True)])
            line10 = LineString([roadline.interpolate(0.45, normalized = True), roadline.interpolate(0.5, normalized = True)])
            line11 = LineString([roadline.interpolate(0.5, normalized = True), roadline.interpolate(0.55, normalized = True)])
            line12 = LineString([roadline.interpolate(0.55, normalized = True), roadline.interpolate(0.6, normalized = True)])
            mid9 = line9.interpolate(0.5, normalized = True)
            mid10 = line10.interpolate(0.5, normalized = True)
            mid11 = line11.interpolate(0.5, normalized = True)
            mid12 = line12.interpolate(0.5, normalized = True)
            #calculate the distance between end point1 of a mid line and the corresponding point in the road line and 
            # the distance between the other end point and the same point in road line
            #if distance differ a lot append them to a list so we can check the len of the list after
            differences = []
            dif8start = LineString([end8one, mid9]).length
            dif8end = LineString([end8tow, mid9]).length
            if (dif8start-dif8end) > 0:
                if (dif8start - dif8end) > 0.7*dif8end:
                    differences.append(1)
            elif (dif8start-dif8end) < 0:
                if dif8end-dif8start > 0.7*dif8start:
                    differences.append(1)

            dif9start = LineString([end9one, mid10]).length
            dif9end = LineString([end9tow, mid10]).length 
            if (dif9start-dif9end) > 0:
                if (dif9start - dif9end) > 0.81*dif9end:
                    differences.append(1)
            elif (dif9start-dif9end) < 0:
                if dif9end-dif9start > 0.81*dif9start:
                    differences.append(1)
            
            dif10start = LineString([end10one, mid11]).length 
            dif10end = LineString([end10tow, mid11]).length
            if (dif10start-dif10end) > 0:
                if (dif10start - dif10end) > 0.81*dif10end:
                    differences.append(1)
            elif (dif10start-dif10end) < 0:
                if dif10end-dif10start > 0.81*dif10start:
                    differences.append(1)
            
            dif11start = LineString([end11one, mid12]).length 
            dif11end = LineString([end11tow, mid12]).length
            if (dif11start-dif11end) > 0:
                if (dif11start - dif11end) > 0.7*dif11end:
                    differences.append(1)
            elif (dif11start-dif11end) < 0:
                if dif11end-dif11start > 0.7*dif11start:
                    differences.append(1)

            if len(differences) >= 2:
                #print("there is a T intersection case2")
                return True
            else:
                #print("the road is wider in the middle but is not T intersection or cross intersection")
                return False
    
        elif size == 1:
            end4one = Point(twentylines[4][1].coords[0])
            end4tow = Point(twentylines[4][1].coords[-1])
            end5one = Point(twentylines[5][1].coords[0])
            end5tow = Point(twentylines[5][1].coords[-1])
            end6one = Point(twentylines[6][1].coords[0])
            end6tow = Point(twentylines[6][1].coords[-1])
            #find the midpoints of the road lines that are in line with the 3 mid lines
            line5 = LineString([roadline.interpolate(0.3, normalized = True), roadline.interpolate(0.4, normalized = True)])
            line6 = LineString([roadline.interpolate(0.4, normalized = True), roadline.interpolate(0.5, normalized = True)])
            line7 = LineString([roadline.interpolate(0.5, normalized = True), roadline.interpolate(0.6, normalized = True)])
            mid5 = line5.interpolate(0.5, normalized = True)
            mid6 = line6.interpolate(0.5, normalized = True)
            mid7 = line7.interpolate(0.5, normalized = True)

            #calculate the distance between end point1 of a mid line and the corresponding point in the road line and 
            # the distance between the other end point and the same point in road line
            #if distance differ a lot append them to a list so we can check the len of the list after
            differences = []
            dif4start = LineString([end4one, mid5]).length 
            dif4end = LineString([end4tow, mid5]).length
            if (dif4start-dif4end) > 0:
                if (dif4start - dif4end) > 0.7*dif4end:
                    differences.append(1)
            elif (dif4start-dif4end) < 0:
                if dif4end-dif4start > 0.7*dif4start:
                    differences.append(1)

            dif5start = LineString([end5one, mid6]).length
            dif5end = LineString([end5tow, mid6]).length
            if (dif5start-dif5end) > 0:
                if (dif5start - dif5end) > 0.81*dif5end:
                    differences.append(1)
            elif (dif5start-dif5end) < 0:
                if dif5end-dif5start > 0.81*dif5start:
                    differences.append(1)
            
            dif6start = LineString([end6one, mid7]).length
            dif6end = LineString([end6tow, mid7]).length 
            if (dif6start-dif6end) > 0:
                if (dif6start - dif6end) > 0.81*dif6end:
                    differences.append(1)
            elif (dif6start-dif6end) < 0:
                if dif6end-dif6start > 0.81*dif6start:
                    differences.append(1)
            
            if len(differences) >= 2:
                print("there is a T intersection case2")
                return True
            else:
                print("the road is wider in the middle but is not T intersection or cross intersection")
                return False
    else :
        return False

def checkforTmaincase(twentylines, roadline, size):
    #print("it passed the test for T intersection of case 2, checking for T intersection Main case!")
    if size == 2: 
        if len(twentylines) > 19:
            first = [twentylines[0][0], twentylines[1][0], twentylines[2][0]]#, twentylines[3], twentylines[4], twentylines[5]]
            last = [twentylines[17][0], twentylines[18][0], twentylines[19][0]]
        elif len(twentylines) > 18:
            first = [twentylines[0][0], twentylines[1][0], twentylines[2][0]] #, twentylines[3], twentylines[4], twentylines[5]]
            last = [twentylines[16][0], twentylines[17][0], twentylines[18][0]]
        elif len(twentylines) > 17:
            first = [twentylines[0][0], twentylines[1][0], twentylines[2][0]] #, twentylines[3], twentylines[4], twentylines[5]]
            last = [twentylines[15][0], twentylines[16][0], twentylines[17][0]]
        else:
            return False
    elif size == 1:
        if len(twentylines) > 9:
            first = [twentylines[0][0], twentylines[1][0], twentylines[2][0]]
            last = [twentylines[8][0], twentylines[9][0]]
        elif len(twentylines)> 8:
            first = [twentylines[0][0], twentylines[1][0]]
            last = [twentylines[7][0], twentylines[8][0]]
        elif len(twentylines)> 7:
            first = [twentylines[0][0], twentylines[1][0]]
            last = [twentylines[6][0], twentylines[7][0]]
        else:
            return False

    f = avgwidth2(first)
    l = avgwidth2(last)
    #print("first lines mean: ", f, " " "and last lines mean: ", l)
    if  f > 2.1*l or l> 2.1*f:
        #print("Main case of T intersection")
        return True
    else:
        #print("it passes the test for T main case intersection")
        return False

def findcrossandx(allpolygons, allpoints4, merg_lines_gdf , epgs):
    print("start to find Cross and X polygons!!")
    plss =list_to_dgf(allpolygons,  epgs)
    pntss = list_to_dgf(allpoints4,  epgs)
    # Among the list with all the polygons, find those polygons that contain at least 1 node of degree 4
    polygons_deg4 = []
    for pp in range(len(plss)):
        polgn = plss.loc[pp, 'geometry']
        p_noholes = close_holes(polgn)
        #if polgn != p_noholes:
            #print("i closed that polygon:")
            #print(polgn)
            #print(p_noholes) 
        neighbor_points = pntss.sindex.query(p_noholes)
        for nei in neighbor_points:
            recc = pntss.loc[nei, 'geometry']
            if p_noholes.contains(recc):
                polygons_deg4.append(polgn)

    print("i have the polygons that contain nodes of degree 4 or higher")
    for p in allpolygons:
        for pp in polygons_deg4:
            if type(p) == list:
                continue
            else:
                if p == pp and p in allpolygons: 
                    allpolygons.remove(p)
                    edges_ofp = merg_lines_gdf.sindex.query(p, predicate= 'intersects')
                    #print(p)
                    #print(len(edges_ofp))
                    if len(edges_ofp)>= 4:
                        #print("weird intersectoin")
                        #print(p)
                        n_p = close_holes(p)
                        allpolygons.append([n_p,"another"])
                    else:
                        #print("cross or x")
                        #print(p)
                        neigss = []
                        for ninn in edges_ofp:
                            neigss.append(merg_lines_gdf.loc[ninn, 'geometry'])
                            #print(merg_lines_gdf.loc[ninn, 'geometry'])
                        if len(neigss) == 2:
                            inter_point = neigss[0].intersection(neigss[1])
                            #print(inter_point)
                            if inter_point.geom_type == 'MultiPoint':
                                #print("weird intersectoin")
                                #print(p)
                                n_p = close_holes(p)
                                allpolygons.append([n_p,"another"])
                            else:
                                l1 = neigss[0].intersection(p)
                                if l1.geom_type == 'MultiLineString' or l1.geom_type == 'GeometryCollection':
                                    #print("the first intersection is multi line ")
                                    #print(l1)
                                    l1 = l1[0]
                                    #print(l1)
                                p1 = Point(l1.coords[0][0], l1.coords[0][1])
                                new_line1 = LineString([inter_point, p1])
                                #print(new_line1)
                                l2 = neigss[1].intersection(p)
                                if l2.geom_type == 'MultiLineString' or l2.geom_type == 'GeometryCollection':
                                    #print("the second intersection is multi line ")
                                    #print(l2)
                                    l2 = l2[0]
                                    #print(l2)
                                p2 = Point(l2.coords[0][0], l2.coords[0][1])
                                new_line2 = LineString([inter_point, p2])
                                #print(new_line2)
                                #check angle between 2 lines
                                x1 = new_line1.coords[1][0] - new_line1.coords[0][0]
                                y1 = new_line1.coords[1][1] - new_line1.coords[0][1]
                                x2 = new_line2.coords[1][0] - new_line2.coords[0][0]
                                y2 = new_line2.coords[1][1] - new_line2.coords[0][1]
                                angle = anglebetvectors(x1,y1,x2,y2)
                                if angle == "w":
                                    #print("cross")
                                    n_p = close_holes(p)
                                    allpolygons.append([n_p,"cross"])
                                else:    
                                    #print(angle)
                                    if angle> 1.4 and angle <1.7:
                                        #print("cross")
                                        n_p = close_holes(p)
                                        allpolygons.append([n_p,"cross"])
                                    else:
                                        #print("X")
                                        n_p = close_holes(p)
                                        allpolygons.append([n_p,"X"])
                        elif len(neigss) == 3:
                            #print("its 3 lines")
                            #print(p)
                            n_p = close_holes(p)
                            allpolygons.append([n_p,"cross"])
                            
                        else:
                            #print(p)
                            n_p = close_holes(p)
                            allpolygons.append([n_p,"another"])
                    

    print("i have cross , x, another_type polygons")
    return allpolygons, plss, pntss

def findtintersections(polygons, int_nodes, epgs, edges, dist,check, another):
    # find all the intersection polygons that are not identified as any type yet (x, cross etc.)
    int_pols = []
    removed = []
    for i in range(len(int_nodes)):
        rec = int_nodes.loc[i,'geometry']
        n = polygons.sindex.query(rec,predicate= "intersects")
        if len(n) == 0:
            continue
        elif len(n) == 1:
            int_p = polygons.loc[n[0], 'geometry']
            #print(int_p)
            if int_p not in removed:
                check.remove(int_p)
                removed.append(int_p)
            int_pols.append(int_p)

    inters = list_to_dgf(int_pols, epgs)
    #for each int polygon, find an edge that crosses it
    #fit the edge to each polygon (fitline), if more than 1 edge crosses the polygon the choose the bigger line
    tss = []
    ass= []
    for intt in range(len(inters)):
        alll = []
        pl1 = inters.loc[intt, 'geometry']
        pl = close_holes(pl1)
        size = 2
        #print(pl)
        fitline_result = fitlinetopolygons(pl, edges)
        fitline = fitline_result[0]
        another_int = fitline_result[1]
        #print(another_int)
        if another_int != 0:
            ass.append(another_int)
        else:
            #print(fitline)
            cutlines = twentlines(fitline)
            for i in cutlines:
                measure_line = measuringline(i, dist, pl, fitline)
                #print(measure_line[1])
                alll.append(measure_line)
            #print(alll)
            if checkforTcase2(alll, fitline , size):
                #print("t")
                tss.append(pl)
                #print(pl)
            elif checkforTmaincase(alll, fitline, size):
                #print("t")
                #print(pl)
                tss.append(pl)
            else:
                #print("a")
                #print(pl)
                ass.append(pl)

    for plgnn in ass:
        another.append(plgnn)
    
    return check, tss, another


#help functions

def avgwidth2(lst):
    s = 0
    for i in lst:
        s += i
    avg = s / len(lst)
    return round(avg, 6)

def makeoffsets(line, dist):
    leftline = line[0].parallel_offset(dist, 'left')
    rightline = line[0].parallel_offset(dist, 'right')
    #print("offsetsss")
    #print(leftline)
    #print(rightline)
    #we find the 2 midpoints of the 2 offset lines in order to connect them and create the prependicular line segment
    midpointleft = leftline.interpolate(0.5, normalized = True)
    midpointright = rightline.interpolate(0.5, normalized = True)
    #print(midpointleft)
    #print(midpointright)
    prependicular_line = LineString([midpointleft, midpointright])
    #print("prep")
    #print(prependicular_line)
    prep_line = [prependicular_line, line[1]]
    return prep_line

def widt(prepedicularline, polygon, roadline):
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
    
    return [w, ln]

def measuringline(line, dist, polygon, fitline):
    prep_line = makeoffsets(line, dist)
    oww = widt(prep_line[0], polygon, fitline)
    width = oww[0]
    if oww[1].geom_type == 'GeometryCollection':
        for mln in oww[1]:
            if mln.geom_type == 'LineString':
                measuringline = mln
                break
    else:
        measuringline = oww[1]
    return [width,measuringline]

def twentlines(line):
    tenlines = []
    line1 = [LineString([line.interpolate(0, normalized = True), line.interpolate(0.05, normalized = True)]), 1]
    line2 = [LineString([line.interpolate(0.05, normalized = True), line.interpolate(0.1, normalized = True)]), 2]
    line3 = [LineString([line.interpolate(0.1, normalized = True), line.interpolate(0.15, normalized = True)]), 3]
    line4 = [LineString([line.interpolate(0.15, normalized = True), line.interpolate(0.2, normalized = True)]), 4]
    line5 = [LineString([line.interpolate(0.2, normalized = True), line.interpolate(0.25, normalized = True)]), 5]
    line6 = [LineString([line.interpolate(0.25, normalized = True), line.interpolate(0.3, normalized = True)]), 6]
    line7 = [LineString([line.interpolate(0.3, normalized = True), line.interpolate(0.35, normalized = True)]), 7]
    line8 = [LineString([line.interpolate(0.35, normalized = True), line.interpolate(0.4, normalized = True)]), 8]
    line9 = [LineString([line.interpolate(0.4, normalized = True), line.interpolate(0.45, normalized = True)]), 9]
    line10 = [LineString([line.interpolate(0.45, normalized = True), line.interpolate(0.5, normalized = True)]), 10]
    line11 = [LineString([line.interpolate(0.5, normalized = True), line.interpolate(0.55, normalized = True)]), 11]
    line12 = [LineString([line.interpolate(0.55, normalized = True), line.interpolate(0.6, normalized = True)]), 12]
    line13 = [LineString([line.interpolate(0.6, normalized = True), line.interpolate(0.65, normalized = True)]), 13]
    line14 = [LineString([line.interpolate(0.65, normalized = True), line.interpolate(0.7, normalized = True)]), 14]
    line15 = [LineString([line.interpolate(0.7, normalized = True), line.interpolate(0.75, normalized = True)]), 15]
    line16 = [LineString([line.interpolate(0.75, normalized = True), line.interpolate(0.8, normalized = True)]), 16]
    line17 = [LineString([line.interpolate(0.8, normalized = True), line.interpolate(0.85, normalized = True)]), 17]
    line18 = [LineString([line.interpolate(0.85, normalized = True), line.interpolate(0.9, normalized = True)]), 18]
    line19 = [LineString([line.interpolate(0.9, normalized = True), line.interpolate(0.95, normalized = True)]), 19]
    line20 = [LineString([line.interpolate(0.95, normalized = True), line.interpolate(1, normalized = True)]), 20]
    tenlines.append(line1)
    tenlines.append(line2)
    tenlines.append(line3)
    tenlines.append(line4)
    tenlines.append(line5)
    tenlines.append(line6)
    tenlines.append(line7)
    tenlines.append(line8)
    tenlines.append(line9)
    tenlines.append(line10)
    tenlines.append(line11)
    tenlines.append(line12)
    tenlines.append(line13)
    tenlines.append(line14)
    tenlines.append(line15)
    tenlines.append(line16)
    tenlines.append(line17)
    tenlines.append(line18)
    tenlines.append(line19)
    tenlines.append(line20)
    return tenlines

def fitlinetopolygons(pl, edges):   
    lines_polygons_pairs = []
    n = edges.sindex.query(pl, predicate = "intersects")
    #print(len(n))
    if len(n) == 2:
        ln1 = edges.loc[n[0], 'geometry']
        ln2 = edges.loc[n[1], 'geometry']
        fitline111 = ln1.intersection(pl)
        fitline222 = ln2.intersection(pl)
        in1111 = fitline111.intersection(fitline222)
        if in1111.geom_type != 'Point':
            another_int = pl
        else:
            another_int = 0
    else:
        another_int = 0
        
    if len(n) > 1:
        lengthss = []
        for nn in n:
            ln = edges.loc[nn, 'geometry']
            inn = ln.intersection(pl)
            if inn.geom_type == 'MultiLineString':
                lng1 = 0
                for kk in inn:
                    lng = kk.length
                    if lng > lng1:
                        lng1 = lng
                        kkk= kk
                lengthss.append([lng1, kkk])
            elif inn.geom_type == 'LineString':
                lengthss.append([inn.length, inn])
            else:
                print("what?")
        
        a = sorted(lengthss, key=lambda x: x[0])
        fitline = a[-1][1]
        #print(fitline)

    else:
        ln = edges.loc[n[0], 'geometry']
        inn = ln.intersection(pl)
        if inn.geom_type == 'MultiLineString':
            lng1 = 0
            for kk in inn:
                lng = kk.length
                if lng > lng1:
                    lng1 = lng 
                    kkk = kk
            fitline = kkk
            #print(fitline)
        elif inn.geom_type == 'LineString':
            fitline = inn
            #print(fitline)      
        else:
            print("what? what?")

    return fitline, another_int

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

def generatestats(p,pntgdf, nm):
    '''Generate some basic stats for a road polygons shapefile:
    Total number of polygong, total_intersection_polygons, percentage of int_polygons, Mean_area (m2), Std from mean, 
    Median_area (m2),Max road polygon area, min road polygons area
    '''
    list_areas = []
    smalls = []
    bigs = []
    for i in range(len(p)):
        pl = p.loc[i, 'geometry']
        list_areas.append(pl.area)

    a = areameanmedian(list_areas)
    mean = a[0]
    median = a[1]
    std = a[2]
    num_pols = len(p)
    max_p = max(list_areas)
    min_p = min(list_areas)

    for i in range(len(p)):
        pl = p.loc[i, 'geometry']
        if pl.area< 15:
            smalls.append(pl)
        
        if pl.area>2000:
            bigs.append(pl)

    percent_smalls = (len(smalls)*100) / len(p)
    percent_bigs = (len(bigs)*100) / len(p)

    ll=[]
    for i in range(len(p)):
        pl = p.loc[i,'geometry']
        n = pntgdf.sindex.query(pl, predicate = 'contains')   
        if len(n)>0:
            ll.append(pl)

    num_int_polygons = len(ll)
    perc_int_pols = (num_int_polygons*100) / len(p)

    with open(nm, 'w', newline='') as fileeee:
        fieldnames = ['total_polygons', 'total_intersection_polygons', '%int_polygons', 'Mean_area (m2)', 'Std', 'Median_area (m2)','Max road polygon area',
         'Min road polygon area', 'Total small polygons', '%small polygons', 'Total big polygons', '%big polygons' ]
        writer = csv.DictWriter(fileeee, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'total_polygons': len(p), 'total_intersection_polygons': len(ll), '%int_polygons':perc_int_pols, 'Mean_area (m2)': mean,
        'Std': std, 'Median_area (m2)':median ,'Max road polygon area': max_p, 'Min road polygon area':min_p, 
        'Total small polygons':len(smalls), '%small polygons':percent_smalls, 'Total big polygons':len(bigs), '%big polygons':percent_bigs})

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
    #print(mean -std)
    #print(mean+std)
    for i in range(len(l)):
        ln = l.loc[i, 'geometry']
        if ln.length < mean - (mean*0.9):
            smalls.append(ln)
        if ln.length > 2*mean:
            bigs.append(ln)
        
    #print(len(smalls))
    #print(len(bigs))
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

def anglebetvectors(x1,y1,x2,y2):
    onom = (x1*x2 + y1*y2)
    paron = (math.sqrt((x1**2 + y1**2)) * math.sqrt((x2**2 + y2**2)))
    if paron == 0:
        return "w"
    else:
        this = onom / paron
        angle =  math.acos(this)
        return angle

def delete_files_2():
    os.remove("projected_lines_shapefile.shp")
    os.remove("projected_lines_shapefile.cpg")
    os.remove("projected_lines_shapefile.dbf")
    os.remove("projected_lines_shapefile.prj")
    os.remove("projected_lines_shapefile.shx")
    os.remove("projected_polygons_shapefile.shp")
    os.remove("projected_polygons_shapefile.cpg")
    os.remove("projected_polygons_shapefile.dbf")
    os.remove("projected_polygons_shapefile.prj")
    os.remove("projected_polygons_shapefile.shx")

def readpolygons(gpdframe):
    allpolygons = []
    for i in range(len(gpdframe)):
        p = gpdframe.loc[i,'geometry']
        validit = p.is_valid
        if validit == False:
            #print(explain_validity(p))
            clean_p = p.buffer(0)
            #print(explain_validity(clean_p))
            allpolygons.append(clean_p)
        elif validit == True:
            allpolygons.append(p)
    return allpolygons

def readlines(shapefile_name):
    polygons_test = shapefile.Reader(shapefile_name)
    shapes = polygons_test.shapes()
    alllines = []
    road_id = 0
    for i in range(len(shapes)):
        shape = shapes[i].points
        p = LineString (shape)
        alllines.append([p, road_id])
        road_id+=1
    return alllines

def readnodes(p_nodes):
    al = gpd.read_file(p_nodes)
    allnodes = []
    for i in range(len(al)):
        allnodes.append(al.loc[i, 'geometry'])
    return allnodes

# PRE DEFINED VALUES
epgs = "EPSG:3879"
bufer = 19.5
dist_for_split = 27
area_for_mrg = 32
dist_for_offsets = 70
polygons_file =  "helsinki_tested_area.shp"
nodes_path = "C:/Users/chcha/Desktop/GEOMATICS/2nd Year/Thesis/p2_technical/code/methodology1/graph/nodes.shp"
edges_path =  "C:/Users/chcha/Desktop/GEOMATICS/2nd Year/Thesis/p2_technical/code/methodology1/graph/edges.shp"

#re-create road polygons 
new_roads = recreate_int_polygons(polygons_file, epgs , bufer, dist_for_split, area_for_mrg, nodes_path, edges_path)

print("POLYGONS HAS BEEN RE-CREATED!")
print("START INTERSECTION IDENTIFICATION!")
print("PRE-PROCESSING!")

#recreated polygons = allpolygons
#check if in the polygons file there are roundabouts and read the final re-created polygons
if new_roads[4] != 'no roundabouts': 
    recreated_polygons = gpd.read_file(new_roads[4])
    allpolygons = readpolygons(recreated_polygons)
else:
    recreated_polygons = gpd.read_file(new_roads[0])
    allpolygons = readpolygons(recreated_polygons)

#Generate stats for the polygons shapefile
int_nodes_all = new_roads[7]
print("GENERATE STATS FOR THE RE-CREATED POLYGONS")
generatestats(recreated_polygons, int_nodes_all,'polygon_stats.csv')

#nodes
p_nodes4 = new_roads[1]
allnodes = new_roads[8]
allpoints4 = readnodes(p_nodes4)
mtr_nodes = new_roads[6]
# these are all the points of degree 4 or higher: allpoints4
# these are all the intersection nodes: allnodes

#edges
edges_shp = new_roads[2]
e = gpd.read_file(edges_shp)
e2 = e.to_crs(epgs)
print("GENERATE STATS FOR CENTERLINES")
generatestats2(e2, 'unmerged_lines_stats.csv')
#merging original centerlines that have the same attribute 'name' and 'highway', 
# this merged lines will just used for several steps (identify intersections) but they wont be linked with width estimations or other information
m_lines = merge_lines_namehighway(edges_shp, epgs)
merg_lines = m_lines[0]
merg_lines_gdf = m_lines[1]
lnsgdf = gpd.read_file(merg_lines)
generatestats2(lnsgdf, 'merged_lines_stats.csv')

alllines = []
cnt = 0
for linee in range(len(lnsgdf)):
    l = [lnsgdf.loc[linee,'geometry'], cnt]
    cnt+=1
    alllines.append(l)
    
#exclude the road polygons that only crossed by motorways or trunks
#Use attribute of the edges (highway)
for popp in allpolygons:
    n = e2.sindex.query(popp, predicate= "intersects")
    tt = []
    if len(n)>= 1:
        for nn in n:
            #print(e.loc[nn,'geometry'])
            rec = e2.loc[nn]
            if rec['highway'] == 'motorway_link' or rec['highway'] == 'motorway' or rec['highway'] == 'trunk_link' or rec['highway'] == 'trunk':
                tt.append(rec)
    if len(tt) != 0 and len(tt) == len(n):
        #print("this polygon is crossed only by motorways!! i remove it from the final width estimation process:")
        #print(popp)
        allpolygons.remove(popp)

#Use attribute of nodes
#mtr_nodes returned from function recreate_int_polygons
if mtr_nodes.empty:
    print("no motor nodes")
else:
    for pngl in allpolygons:
        n = e2.sindex.query(pngl)
        n = mtr_nodes.sindex.query(pngl)
        if len(n) >= 1:
           # print("this polygon contains a motor way node, i remove it: ")
            #print(pngl)
            allpolygons.remove(pngl)

print("these are the final polygons (before intersection identification):  ")
print(len(allpolygons))          
polss = list_to_dgf(allpolygons, epgs)


#find x and cross intersection polygons
tss = findcrossandx(allpolygons, allpoints4, merg_lines_gdf, epgs)
allpolygons = tss[0]
plss = tss[1]
pntss = tss[2]

print("start find roundabout polygons")
roundaboutss = new_roads[5]
if roundaboutss != "no roundabouts":
    for r in roundaboutss:
        nn = plss.sindex.query(r)
        for nin in nn:
            p = plss.loc[nin, 'geometry']
            if r == p and p in allpolygons:
                allpolygons.remove(p)
                allpolygons.append([p,"round"])

print("i have roundabout polygons")
cross = []
rounds = []
x = []
another = []  
check = []
for iiii in allpolygons:
    if type(iiii) == list:
        if iiii[1] == "round":
            rounds.append(iiii[0])
        elif iiii[1] == "cross":
            cross.append(iiii[0])
        elif iiii[1] == "X":
            x.append(iiii[0])
        elif iiii[1] == "another":
            another.append(iiii[0])
    else:
        #it is either road edge polygon or intersection of another type (T, double T, Y)
        check.append(iiii)

ah = list_to_dgf(check, epgs)
ah2 = list_to_dgf(allnodes,epgs)

final_test = findtintersections(ah, ah2, epgs,lnsgdf,dist_for_offsets,check, another)

final_round_edges = final_test[0]
another = final_test[2]
t_ints = final_test[1]

df_roads = list_to_dgf(final_round_edges, epgs)
df_roads.to_file("standardized_road_edges.shp")

df_t = list_to_dgf(t_ints, epgs)
df_t.to_file('T_polygons.shp')

df_cross = list_to_dgf(cross,epgs)
df_cross.to_file("cross_polygons.shp")

df_x = list_to_dgf(x,epgs)
df_x.to_file("X_polygons.shp")

df_rounds = list_to_dgf(rounds,epgs)
df_rounds.to_file("round_polygons.shp")

df_another = list_to_dgf(another,epgs)
df_another.to_file("another_type_polygons.shp")

