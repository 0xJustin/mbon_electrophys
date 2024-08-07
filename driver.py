import neuprint_functions
import numpy as np
import json
import pandas

def main(neuron_id, supertype='KC'):
    types, bodyids = neuprint_functions.get_all_neurons_presynaptic_to(neuron_id, supertype=supertype)
    synapses = neuprint_functions.pull_mb_synapses(bodyids, post_id=neuron_id, only_mb=True, n_iter=4)
    # Factor to adjust from pixels to um coordinates
    coordinates, lengths = neuprint_functions.pre_syn_coord_dict(synapses, bodyids)
    print("The Total number of PRE>POST synapses is = "+str(len(synapses))) 
    print("The total number of presynaptic KCs is = "+str(len(coordinates))) 
    print("The average number of KC>MBON synapses per KC is = "+str(np.mean(lengths)))
    coordinates_json = json.dumps(coordinates, indent = 4)
    with open(f'synapse_coordinates_scaled_{neuron_id}.json', 'w') as output:
        output.write(coordinates_json)
    output.close()

    skel_graph = neuprint_functions.pull_skeleton(neuron_id, skel_format='nx')
    skel = neuprint_functions.pull_skeleton(neuron_id, skel_format='pandas')
    skel_swc, relabel_dict = neuprint_functions.relabel_skeleton_swc(skel_graph, skel)
    with open("neuron_skeleton.swc", "w") as output:
        output.write(skel_swc)
    return skel_swc
