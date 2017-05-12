import argparse
import csv
import datetime
import numpy as np
import itertools
import time
import math
from collections import Counter
from operator import itemgetter

LAT_OFFSET = 40.5000000
LONG_OFFSET = -74.25000000
DATE_OFFSET = datetime.date(2015,1,1).toordinal()

def weight_tuple(raw_tuple, cell_size, timestep):
    DEBUG = False
    x = int((raw_tuple[0] - LAT_OFFSET) / cell_size)
    y = int((raw_tuple[1] - LONG_OFFSET) / cell_size)
    # probably too much done in one step...
    raw_date = raw_tuple[2].split(" ")[0].split("-")
    # probably too much done in one step...
    z = int((datetime.date(int(raw_date[0]),int(raw_date[1]),int(raw_date[2])).toordinal() - DATE_OFFSET) / timestep)
    if DEBUG:
        print((x,y,z))
    return (x,y,z)

def calculate_G_score(cell_values, number_of_cells):
    DEBUG = False
    n = number_of_cells
    ## calculate mean
    mean = x_bar = sum(cell_values) / n
    ## calculate S
    S = math.sqrt(((sum([ x * x for x in cell_values])) / n ) - (math.pow(mean,  2)))

    numerator  = sum(cell_values) - (mean * 27.0) # I think the Wij is const, ie. 1
    denominator = S * math.sqrt(((n * 27.0) - math.pow(27, 2)) / (n - 1))
    if denominator == 0.0:
        return 0.0 #This isn't failure proof, it's just to prevent crash
    if DEBUG:
        print((n, mean, S, numerator, denominator))
    return numerator / denominator

def unweight_tuple(raw_tuple, cell_size, timestep):
    x = raw_tuple[0] * cell_size + LAT_OFFSET
    y = raw_tuple[1] * cell_size + LONG_OFFSET
    z = raw_tuple[2]
    return (x,y,z)

def main():
    DEBUG = False
    parser = argparse.ArgumentParser()
    parser.add_argument('taxi_data_file')
    parser.add_argument('output_file')
    parser.add_argument('cell_size')
    parser.add_argument('timestep')

    print(time.time())

    args = parser.parse_args()
    cell_size = float(args.cell_size)
    timestep = int(args.timestep)
    point_list = []

    with open(args.taxi_data_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            #turn the row into a (lat,long,time) tuple
            if row['dropoff_latitude']:
                x = float(row['dropoff_latitude'])
            if row['dropoff_longitude']:
                y = float(row['dropoff_longitude'])
            if row['tpep_dropoff_datetime']:
                z = row['tpep_dropoff_datetime']
            if x != 0.0 and y != 0.0:
                weighted_tuple = weight_tuple((x,y,z),cell_size,timestep)
                point_list.append(weighted_tuple)


    # Consider cleaning off the junk coordinates (in the water, outside the 5 burroughs envelope)

    # Count the number of events happening at each coordinate Lat,Long,Time
    # consider breaking this into a func
    histogram = dict(Counter(point_list))

    # intiialize the 3D time cube matrix
    # consider breaking this into a func
    max_lat = max(histogram.keys(), key=lambda x: x[0])[0]
    max_long = max(histogram.keys(), key=lambda x: x[1])[1]
    max_timestep = max(histogram.keys(), key=lambda x: x[2])[2]

    if DEBUG:
        print((max_lat,max_long,max_timestep))

    # loop through the whole matrix
    # calculate the Gscore of each point by looking at its neighbors
    cube = np.array(list(itertools.product((-1,0,1), (-1,0,1), (-1,0,1))))

    G_score_out = []
    cell_count = max_lat * max_long * max_timestep

    for x,y,z in itertools.product(range(max_lat), range(max_long), range(max_timestep)):
        if DEBUG:
            print((x,y,z))
            print(cube + (x,y,z))
        cell_vals = []


        cell_and_neighbors = tuple([tuple(row) for row in cube + (x,y,z)])
        if DEBUG:
            print(cell_and_neighbors)
        for cell in cell_and_neighbors:
            cell_vals.append(histogram.get(cell,0))

        G_score_out.append([(x,y,z),calculate_G_score(cell_vals,cell_count)])

    G_score_out.sort(key=itemgetter(1),reverse=True)
    top_50 = G_score_out[:50]

    top_50_out = []

    for score in top_50:
        top_50_out.append([unweight_tuple(score[0], cell_size, timestep), score[1]])

    with open(args.output_file, 'w+', newline='') as outfile:
        writer = csv.writer(outfile)
        #this may take some additional post processing... =/
        writer.writerows(top_50_out)
    print(time.time())

if __name__ == '__main__':
    main()
