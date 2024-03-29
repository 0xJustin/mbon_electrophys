import neuprint
import json
import pandas as pd
from neuprint import Client
import networkx as nx
import numpy as np

token = json.load(open("neuprint_token.json"))["token"]
if token == "":
    raise Exception("You need to add your neuprint token to neuprint_token.json")
try:
    c = Client('neuprint.janelia.org', dataset='hemibrain:v1.2.1', token=token)
except:
    print("Failed to connect to neuprint. Check your token and make sure you are connected to the internet.")

def get_types_and_bodyIds(supertype="KC"):
    """
    Returns two lists: a list of neuron types and a list of body IDs for neurons of the given supertype.
    :param client: A neuprint client object.
    :param supertype: A string representing the supertype of neurons to query. Default is "KC".
    :return: A tuple of two lists: a list of neuron types and a list of body IDs.
    """

    # This is called a cypher query. It's a way of querying the neuprint database precisely
    # https://neuprint.janelia.org/help/cypherexamples?dataset=hemibrain%3Av1.2.1 go here for more examples
    where = f"""(a.type CONTAINS "{supertype}")"""
    q = f" MATCH (a :`hemibrain_Neuron`) WHERE {where} RETURN a.bodyId, a.type"
    lh_table = neuprint.fetch_custom(q)
    return lh_table["a.type"].tolist(), lh_table["a.bodyId"].tolist()
    
def get_all_neurons_presynaptic_to(bodyId, supertype="KC"):
    """
    Returns two lists: a list of neuron types and a list of body IDs for all neurons presynaptic to the given body ID.
    :param client: A neuprint client object.
    :param bodyId: A string representing the body ID of the postsynaptic neuron.
    :param supertype: A string representing the supertype of neurons to query. Default is "KC".
    :return: A tuple of two lists: a list of neuron types and a list of body IDs.
    """

    # This is called a cypher query. It's a way of querying the neuprint database precisely
    # https://neuprint.janelia.org/help/cypherexamples?dataset=hemibrain%3Av1.2.1 go here for more examples
    where = f"""(a.type CONTAINS "{supertype}" AND b.bodyId = {bodyId})"""
    q = f"""MATCH (a:`hemibrain_Neuron`)-[w:ConnectsTo]->(b:`hemibrain_Neuron`) WHERE {where} RETURN a.bodyId, a.type"""
    lh_table = neuprint.fetch_custom(q)
    return lh_table["a.type"].tolist(), lh_table["a.bodyId"].tolist()

def pull_synapses(pre_ids=None, post_id=612371421, n_iter=4, pre_or_post="pre"):
    """
    Returns a pandas DataFrame containing synapse information for the given presynaptic neurons and postsynaptic neuron.
    :param client: A neuprint client object.
    :param pre_ids: A list of strings representing the body IDs of the presynaptic neurons.
    :param post_id: A string representing the body ID of the postsynaptic neuron. Default is 612371421.
    :param only_mb: A boolean indicating whether to only include synapses within the mushroom body. Default is True.
    :param n_iter: An integer representing the number of iterations to use when querying neuprint. Default is 4.
    :return: A pandas DataFrame containing synapse information.
    """
    # Uncomment this line if you only want synapses who are geographically inside the mushroom body
    synapses_df = pd.DataFrame()

    if pre_or_post == 'both':
        directions = ['pre', 'post']
    else:
        directions = [pre_or_post]

    # Replace with your MBON body ID
    # IMPORTANT - make sure you have the right hemisphere
    # Hemibrain only has segmented the right hemisphere, but PORTIONS
    # of the left hemisphere are included in the dataset, including MBONs
    # You can see this by going to the filters when you are searching for a neuron in neuprint
    # And checking the "Instance" box- you will see that some neurons are "_L" and some are "_R
    n_neurons = len(pre_ids)
    # Iterating in a list because neuprint has a timeout (maximum time its servers will crunch on your request) that is pretty low
    # So it works bettter if you break the request into a series of smaller requests
    for i in range(n_iter):
        # To break up the requests, I just divide up the list of KCs
        if not pre_ids is None:
            presynaptic_neurons = pre_ids[i * n_neurons // n_iter : (i + 1) * n_neurons // n_iter]
        
        for direction in directions:
            if direction == 'pre':
                pre = presynaptic_neurons if not pre_ids is None else None
                post = post_id
            else:
                pre = post_id
                post = presynaptic_neurons if not pre_ids is None else None

            presynaptic_criteria = neuprint.NeuronCriteria(bodyId=pre)
            postsynaptic_criteria = neuprint.NeuronCriteria(bodyId=post)
            try:
                synapses = neuprint.fetch_synapse_connections(presynaptic_criteria, postsynaptic_criteria, client=c)
            except Exception as e:
                print("ERROR: ", e)
                print("Failed to get synapses- probably a neuprint timeout, try increasing n_iter")
                
            synapses_df = synapses_df._append(synapses)
    synapses_df.loc[:, "pre_or_post"] = "post"
    synapses_df.loc[synapses_df["bodyId_pre"] == post_id, "pre_or_post"] = "pre"
    return synapses_df.reset_index()

def pre_syn_coord_dict(synapses, bodyids, factor=8*(10**(-3))):
    coordinates = {}
    for ID in bodyids:
        coordinates[str(ID)] = []
    for i in range(len(synapses)):
        x = int(synapses.loc[i, 'x_post'])
        y = int(synapses.loc[i, 'y_post'])
        z = int(synapses.loc[i, 'z_post'])
        x_scaled = factor * x
        y_scaled = factor * y
        z_scaled = factor * z
        coordinates[str(synapses.loc[i,'bodyId_pre'])].append([x_scaled, y_scaled, z_scaled])

    lengths = []
    for ID, coords in coordinates.items():
        lengths.append(len(coords))
    return coordinates, lengths

def pull_skeleton(body_id, healdist=200, return_kdtree=False, skel_format="pandas"):
    # Skeleton pulling and healing code
    # Another parameter you might like to include is "with_distances=True"
    """
    Fetches the skeleton for the given body ID and returns it as a pandas DataFrame.
    :param client: A neuprint client object.
    :param body_id: A string representing the body ID to fetch the skeleton for.
    :param healdist: An integer representing the distance to use for healing the skeleton. Default is 200.
    :param return_kdtree: A boolean indicating whether to return the skeleton as a cKDTree. Default is False.
    :return: A pandas DataFrame representing the skeleton.
    """
    skel_mbon = neuprint.skeleton.fetch_skeleton(body_id, client=c, heal=healdist, format=skel_format)
    # From https://connectome-neuprint.github.io/neuprint-python/docs/client.html#neuprint.client.Client.fetch_skeleton:
    #  If True, a ‘distance’ column (or edge attribute) will be added to the dataframe (or nx.Graph), indicating the distances from each node to its parent node.

    # You might also find it useful to interpret the skeleton as a cKDTree
    # Basically, KD trees make it much easier to find the points "Close" to a given point
    # For analyzing many many synapse locations, this can dramatically speed up your code
    # https://en.wikipedia.org/wiki/K-d_tree
    if return_kdtree:
        from scipy.spatial import cKDTree
        tree = cKDTree(list(zip(skel_mbon["x"], skel_mbon["y"], skel_mbon["z"])))
        return skel_mbon, tree
    return skel_mbon

def bfs_relabel(skel):
    """
    Performs a breadth-first relabeling on a list of nodes with parent indices.
    :param node_list: A list of nodes, where each node is a tuple of the node index and its parent index.
    :return: A list of node indices visited during the search.
    """
    relabel = {}
    visited = []
    queue = [1]
    counter = 1 
    while queue:
        node_index = queue.pop(0)
        if node_index not in visited:
            visited.append(node_index)
            neighbors = list(skel.predecessors(node_index))
            queue.extend(neighbors)
            relabel[node_index] = counter
            counter += 1
    return relabel

def relabel_skeleton_swc(skel_graph, skel_df, row_name="rowId", link_name="link"):
    relabel_dict = bfs_relabel(skel_graph)
    relabel_dict[-1] = -1
    skel_df[row_name] = skel_df[row_name].map(relabel_dict)
    skel_df[link_name] = skel_df[link_name].map(relabel_dict)
    swc = "# rowId type x y z radius link\n"
    for i in range(skel_df.shape[0]):
        s = f"{int(skel_df.iloc[i][row_name])} 0 {skel_df.iloc[i].x} {skel_df.iloc[i].y} {skel_df.iloc[i].z} {skel_df.iloc[i].radius} {int(skel_df.iloc[i][link_name])}\n"
        swc += s
    return swc, relabel_dict

def match_synapses_to_tree(synapses, skel_graph, relabel_dict):
    from scipy.spatial import cKDTree

    skel_graph = nx.relabel_nodes(skel_graph, relabel_dict)
    skel_coords = np.array([[skel_graph.nodes[i]["x"], skel_graph.nodes[i]["y"], skel_graph.nodes[i]["z"]] for i in skel_graph.nodes])
    synapses_coords = np.array([[synapses.loc[i, 'x_post'], synapses.loc[i, 'y_post'], synapses.loc[i, 'z_post']] if synapses.loc[i, 'pre_or_post'] == 'post' else [synapses.loc[i, 'x_pre'], synapses.loc[i, 'y_pre'], synapses.loc[i, 'z_pre']] for i in synapses.index])
    # Create a cKDTree from synapses_coords
    skel_tree = cKDTree(skel_coords)
    # Calculate the distance matrix using synapses_tree and skel_coords
    dm = skel_tree.query(synapses_coords)[1]
    node_to_synapses = {}
    for i, node in enumerate(skel_graph.nodes()):
        # Find the indices of synapses closest to the current node
        indices = np.where(dm == i)[0]
        # Map the current node to the list of indices
        node_to_synapses[node] = indices

    return dm, node_to_synapses