import neuprint
import json
import pandas as pd
from neuprint import Client

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
    where = f"""(a.type CONTAINS "{supertype}")"""
    q = f" MATCH (a :`hemibrain_Neuron`)-[w :ConnectsTo]->(b :`hemibrain_Neuron`) WHERE {where} AND b.bodyId = \"{bodyId}\" RETURN a.bodyId, a.type"
    lh_table = neuprint.fetch_custom(q)
    return lh_table["a.type"].tolist(), lh_table["a.bodyId"].tolist()

def pull_mb_synapses(pre_ids, post_id=612371421, only_mb=True, n_iter=4):
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
    if only_mb:
        is_MB = """AND NOT apoc.convert.fromJsonMap(w.roiInfo)["MB(+ACA)(R)"] IS NULL"""
    synapses_df = pd.DataFrame()

    # Replace with your MBON body ID
    # IMPORTANT - make sure you have the right hemisphere
    # Hemibrain only has segmented the right hemisphere, but PORTIONS
    # of the left hemisphere are included in the dataset, including MBONs
    # You can see this by going to the filters when you are searching for a neuron in neuprint
    # And checking the "Instance" box- you will see that some neurons are "_L" and some are "_R
    n_neurons = len(pre_ids)
    # Iterating in a list because neuprint has a timeout (maximum time its servers will crunch on your request) that is pretty low
    # So it works bettter if you break the request into a series of smaller requests
    n_iter = 4
    for i in range(n_iter):
        # To break up the requests, I just divide up the list of KCs
        presynaptic_neurons = pre_ids[i * n_neurons // n_iter : (i + 1) * n_neurons // n_iter]
        presynaptic_criteria = neuprint.NeuronCriteria(bodyId=presynaptic_neurons)
        postsynaptic_criteria = neuprint.NeuronCriteria(bodyId=post_id)
        try:
            synapses = neuprint.fetch_synapse_connections(presynaptic_criteria, postsynaptic_criteria, client=client)
        except Exception as e:
            print("ERROR: ", e)
            print("Failed to get synapses- probably a neuprint timeout, try increasing n_iter")
            
        synapses_df = synapses_df._append(synapses)
    return synapses_df

def pull_skeleton(body_id, healdist=200, return_kdtree=False):
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
    skel_mbon = neuprint.skeleton.fetch_skeleton(body_id, client=client, heal=healdist)
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

def graph_from_skeleton(skel):
    """
    Converts a skeleton pandas DataFrame to a networkx graph object.
    :param skel: A pandas DataFrame representing the skeleton.
    :return: A networkx graph object representing the skeleton.
    """
    import networkx as nx
    skeleton = skel.to_undirected()
    return skeleton

def relabel_swc(skel):
    return skel