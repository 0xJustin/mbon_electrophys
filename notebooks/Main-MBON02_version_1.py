# Import python library to access neuprint database
import neuprint
from neuprint import Client, fetch_synapses, fetch_synapse_connections, SynapseCriteria as SC
from neuron import h, gui
import csv
import pandas as pd
import statistics
import json
import numpy as np
import random
import itertools
from scipy.spatial import distance
import matplotlib.pyplot as plt
import statistics
import copy

# Creating your personal client to fetch neuprint dataset with toke
c = Client('neuprint.janelia.org', dataset='hemibrain:v1.2.1',
token='eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6Imtob3Noa2hvdW1haHNhQGdtYWlsLmNvbSIsImxldmVsIjoibm9hdXRoIiwiaW1hZ2UtdXJsIjoiaHR0cHM6Ly9saDMuZ29vZ2xldXNlcmNvbnRlbnQuY29tL2EvQUNnOG9jSkwwemFsaDlTUU1tYTA1ZTVpeEhVWTlXRjNoNkhrUFFzUlZ3QkFwLU1mPXM5Ni1jP3N6PTUwP3N6PTUwIiwiZXhwIjoxODc0ODA1MTQyfQ.hY_wfjSAbChhP6C3C7KpvxLan7_9yh5hRYuAr-uN3nk')
c.fetch_version()
version = c.fetch_version()
print('Version = '+str(version))

# Fetching data from neuprint

where = f"""(a.type CONTAINS "KC")"""
q = f" MATCH (a :`hemibrain_Neuron`) WHERE {where} RETURN a.bodyId, a.type"
lh_table = neuprint.fetch_custom(q)
MBON02_bodyId = 424789697
n_kcs = len(lh_table["a.bodyId"].tolist())
synapses_df = pd.DataFrame()
n_iter = 6
for i in range(n_iter):
    presynaptic_neurons = lh_table["a.bodyId"].tolist()[i * n_kcs // n_iter : (i + 1) * n_kcs // n_iter]
    presynaptic_criteria = neuprint.NeuronCriteria(bodyId=presynaptic_neurons)
    postsynaptic_criteria = neuprint.NeuronCriteria(bodyId=MBON02_bodyId)
    synapses = neuprint.fetch_synapse_connections(presynaptic_criteria, postsynaptic_criteria, client=c)
    synapses_df = synapses_df._append(synapses)
synapses_df.to_csv('KC-MBON02.csv')

coords_x = synapses_df.x_pre.astype(float)

# Heal the neuron
#gradually increase heal param until 'more than one tree', then increase more
swc = c.fetch_skeleton(424789697, heal = 200, export_path = 'healed-MBON-02-structure.swc', format='swc')

### Resort the healed Janelia file to meet the NEURON constraint that parent id< child id
# Load healed swc morphology file
with open('healed-MBON-02-structure.swc','r') as data:
  f = data.readlines()[1:]      # Removing the file header

# Process and load the data into a list "original"
original = []

for line in f:
    line = line.strip()
    columns = line.split()
    columns = [float(x) for x in columns]
    original.append(columns)

print('The length of original is = '+str(len(original)))
print('The entry length of original is '+str(statistics.mean([len(line) for line in original])))

# Now need to expand each line to give a column for new index and new connection partner
# new index | original index | ... | original connection | new connection
extended = copy.deepcopy(original)
for line in extended:
    line.insert(0, line[0])
    line.insert(-1, line[-1])
    line[0] = int(line[0])
    line[1] = int(line[1])
    line[-2] = int(line[-2])
    line[-1] = int(line[-1])

print('The length of extended is = '+str(len(extended)))
print('The entry length of extended is '+str(statistics.mean([len(line) for line in extended])))

sifted = []
flipped_order = []
others = []

n = 1
for line in extended:
    mod_line = copy.deepcopy(line)
    # Check if it's the root
    if line[-2] == -1:
        mod_line[0] = int(len(sifted)+1)
        sifted.append(mod_line)
        n += 1
        continue
    if (line[1] > line[-2]) and (line[-2] in [sifted[i][1] for i in range(0,len(sifted))]):
        mod_line[0] = int(len(sifted)+1)
        index = [i for i, value in enumerate(sifted) if value[1] == line[-2]]
        mod_line[-1] = sifted[index[0]][0]
        sifted.append(mod_line)
        continue
    if (line[1] < line[-2]):
        flipped_order.append(line)
        continue
    else:
        others.append(line)

# Port sorting over to a function
def sort(inputlist, sortout, reversedout, miscout):
    for line in inputlist:
        mod_line = copy.deepcopy(line)
        # Check if it's the root
        if line[-2] == -1:
            mod_line[0] = int(len(sortout)+1)
            sortout.append(mod_line)
            continue
        if (line[1] > line[-2]) and (line[-2] in [sortout[i][1] for i in range(0,len(sortout))]):
            mod_line[0] = int(len(sortout)+1)
            index = [i for i, value in enumerate(sortout) if value[1] == line[-2]]
            mod_line[-1] = sortout[index[0]][0]
            sortout.append(mod_line)
            continue
        if (line[1] < line[-2]):   # if it is flipped in the original, see if connection was already added to sorted list
            if line[-2] in [sortout[i][1] for i in range(0,len(sortout))]:
                mod_line[0] = int(len(sortout)+1)
                index = [i for i, value in enumerate(sortout) if value[1] == line[-2]]
                mod_line[-1] = sortout[index[0]][0]
                sortout.append(mod_line)
                continue
            else:  #if not in sorted list, then save for later
                reversedout.append(line)
                continue
        else:
            miscout.append(line)

sifted = []
flipped_order = []
others = []

sort(extended, sifted, flipped_order, others)

print('The length of sifted is = '+str(len(sifted)))     # 14
print('The length of flipped_order is = '+str(len(flipped_order)))   # 265
print('The length of others is = '+str(len(others)))

# Deal with reversed order ones
flipped_order.reverse()

flipped_order2 = []
others2 = []
sort(flipped_order, sifted, flipped_order2, others2)
print('Second sort of flipped_order (1) ')
print('The length of sifted is = '+str(len(sifted)))
print('The length of flipped_order2 is = '+str(len(flipped_order2)))    #15
print('The length of others2 is = '+str(len(others2)))

flipped_order3 = []
others3 = []
sort(others, sifted, flipped_order3, others3)
print('Third sort of others (1) ')
print('The length of sifted is = '+str(len(sifted)))
print('The length of flipped_order3 is = '+str(len(flipped_order3)))    #
print('The length of others3 is = '+str(len(others3)))

flipped_order4 = []
others4 = []
sort(flipped_order2, sifted, flipped_order4, others4)
print('Fourth sort of flipped_order (2)')
print('The length of sifted is = '+str(len(sifted)))
print('The length of flipped_order4 is = '+str(len(flipped_order4)))
print('The length of others4 is = '+str(len(others4)))

flipped_order5 = []
others5 = []
sort(others3, sifted, flipped_order5, others5)
print('Fifth sort of others (3)')
print('The length of sifted is = '+str(len(sifted)))
print('The length of flipped_order5 is = '+str(len(flipped_order5)))
print('The length of others5 is = '+str(len(others5)))

# So all points are now sifted!!!
print('Check that the new indices are increasing')
intended_indices = [*range(1,34244)]
actual_indices = [i[0] for i in sifted]
if intended_indices == actual_indices:
    print("The indices are all good!")
else:
    print("There is a problem!")

# Check that the new partners are all good. For each point: point(old connection)
# Find the line that has line(old index) = point(old connection)
# Confirm that line(new index) = point(new connection)
good = 0
problems = []
problem_index = []
for point in sifted:
    if point[-1] == -1:
        good += 1
        continue
    if point[-1] > 0:
        partner_index = [i for i, value in enumerate(sifted) if point[-2] == value[1]]
        if point[-1] == sifted[partner_index[0]][0]:
            good += 1
            continue
        else:
            print('The new index and connection dont match at')
            print(good)
            problems.append(point)
            good += 1
            continue
    else:
        print('A point doesnt have a new index that is numeric')
        problem_index.append(point)
        good += 1
        continue

print('If no error messages from connection comparison, look at good number should be 34243')
print(good)

# Trim off the old indices and connection partners....no longer needed
final = copy.deepcopy(sifted)
for line in final:
    del line[1]
    del line[-2]

print('The length of final is = '+str(len(original)))   # Should be 19950
print('The entry length of final is '+str(statistics.mean([len(line) for line in final])))  #SHould be 7

# Output to file
out = open("MBON-02-200-healed-repaired.swc", "w")
for line in final:
    line_string = [str(line[i]) for i in range(len(line))]
    new_line = " ".join(line_string)
    new_line2 = new_line+"\n"
    out.writelines(new_line2)
out.close()

print("Finished writing file")


# adjusting coordinates units
# Take the healed/repaired swc file and convert from pixel coordinates to um
# Each pixel is 8 nm
# pixel x 8*10^(-3) will convert to um

with open('MBON-02-200-healed-repaired.swc','r') as data:
  f = data.readlines()
  # Process and load data into a list "original"
original = []

for line in f:
    line = line.strip()
    columns = line.split()
    columns = [float(x) for x in columns]
    original.append(columns)

print('First 10 lines of file')
for i in range(10):
    print(original[i])
    #print(len(original[i]))

print(int(len(original)))

# Each line has 7 entries defined as:
# line[0] = compartment index
# line[1] = structure type identifier (0 is undefined... all lines in this file)
# line[2] = x coordinate
# line[3] = y coordinate
# line[4] = z coordinate
# line[5] = radius
# line[6] = parent compartment
# So, line 2,3,4,5 all need to be scaled appropriately and lines 0, 1, and 6 must be integers


factor = 8*(10**(-3))
scaled = copy.deepcopy(original)
for line in scaled:
    line[0] = int(line[0])
    line[1] = int(line[1])
    line[2] = factor*line[2]
    line[3] = factor*line[3]
    line[4] = factor*line[4]
    line[5] = factor*line[5]
    line[6] = int(line[6])

print('First 10 lines of scaled file')
for i in range(10):
    print(scaled[i])

print(int(len(scaled)))

# Output to file
out = open('MBON-02-200-Janelia-Scaled.swc','w')
for line in scaled:
    line_string = [str(line[i]) for i in range(len(line))]
    new_line = " ".join(line_string)
    new_line2 = new_line+"\n"
    out.writelines(new_line2)
out.close()

print("Finished writing file")

KC_IDs = list(synapses_df.bodyId_pre.unique())
syn_KCs = fetch_synapse_connections(source_criteria=KC_IDs, target_criteria=MBON02_bodyId)
# Want to collect relevant coordinates in a dictionary sorting by KC
coordinates = {}
for KC in KC_IDs:
    coordinates[str(KC)] = []
# can select a row using loc:
syn_KCs.loc[5,:]
# Factor to adjust from pixels to um coordinates
factor = 8*(10**(-3))

for i in range(len(syn_KCs)):
    x = int(syn_KCs.loc[i, 'x_post'])
    y = int(syn_KCs.loc[i, 'y_post'])
    z = int(syn_KCs.loc[i, 'z_post'])
    x_scaled = factor * x
    y_scaled = factor * y
    z_scaled = factor * z
    coordinates[str(syn_KCs.loc[i,'bodyId_pre'])].append([x_scaled, y_scaled, z_scaled])
lengths = []
for ID, coords in coordinates.items():
    lengths.append(len(coords))

print("The Total number of KC>MBON synapses is = "+str(len(syn_KCs)))   # Should be 12770
print("The total number of presynaptic KCs is = "+str(len(coordinates)))  # Should be 948
print("The average number of KC>MBON synapses per KC is = "+str(statistics.mean(lengths)))   # Should be 13.47

# All looks good, now just need to export synapses as json for use in model
# SAVE THE SYNAPSE DATA
# Put coordinates dictionary into json
coordinates_json = json.dumps(coordinates, indent = 4)
# Export to file
with open('synapse_coordinates_scaled.json', 'w') as output:
    output.write(coordinates_json)
output.close()
print('Finished writing synapses to file')

#saving data to csv
def write_csv(name,data):
  with open(name, 'a') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(data)

def instantiate_swc(filename):
    '''
    Load swc file and instantiate it as cell
    Code source: https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3257
    '''

    # load helper library, included with Neuron
    h.load_file('import3d.hoc')

    # load data
    cell = h.Import3d_SWC_read()
    cell.input(filename)

    # instantiate
    i3d = h.Import3d_GUI(cell,0)
    i3d.instantiate(None)
# Other mbon with DAS (some connectivity issues but doesn't matter for this)
cell = instantiate_swc('MBON-02-200-Janelia-Scaled.swc')

# Model Specification

# Import function automatically creates dendritic sections, named "dend_0[i]" for i from [0, 6143] inclusive
# Create SectionList containing all dendritic sections

dends = h.SectionList()
for sec in h.allsec():
    dends.append(sec=sec)

# Convert SectionList into a Python list dendspy
# Individual sections can be selected by their index i

dendspy = [sec for sec in dends]

# Generationg the data file "synapse-sections.json“
# based on the MBON morphology and the synapse coordinate data. Pairs each synapse
# location with the closest section. Takes a couple of hours to run.

# Collect dictionary of all points defined for each section
MBON_coordinates = {}

for sec in range(len(dendspy)):
    MBON_coordinates[str(sec)] = []
    i = dendspy[sec].n3d()
    for num in range(0, i):
        coord = [dendspy[sec].x3d(num), dendspy[sec].y3d(num), dendspy[sec].z3d(num)]
        MBON_coordinates[str(sec)].append(coord)

# Import synapse coordinates (in um)
with open('synapse_coordinates_scaled.json') as coords:
    synapse_coordinates = json.load(coords)

# For a given point, find the closest neighbor within the list of MBON coordinates
def closest_section(synapse_coordinate):
    dist = None
    sec = None
    for mbon_sec, mbon_coords in MBON_coordinates.items():
        for pt in mbon_coords:
            current_distance = distance.euclidean(pt, synapse_coordinate)
            if dist == None or current_distance < dist:
                dist = current_distance
                sec = mbon_sec
    return sec

print("Finding the closest section to each synapse location")
index = 0  # to keep track of progress of the for loop
synapse_sections = {}  # KC synapses on which sections?

for KC, synapses in synapse_coordinates.items():
    synapse_sections[KC] = []
    index += 1
    print(index)
    for syn_loc in synapses:
        mb_sec = closest_section(syn_loc)
        synapse_sections[KC].append(mb_sec)

print("Now saving the dictionary to a file")
synapses_json = json.dumps(synapse_sections, indent=4)

with open('synapses-sections.json', 'w') as outfile:
    outfile.write(synapses_json)
    outfile.close()

print("Synapses finished saving to the file")


