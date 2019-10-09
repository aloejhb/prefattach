import time
import random
import numpy as np
import networkx as nx
import collections
import matplotlib.pyplot as plt
from scipy.special import gammaln
from operator import mul
from functools import reduce


def introduce_time(x, m):
    return x * 2*m


def total_num_node(t, m):
    return np.floor(t/(2*m)) + 1


def total_adj_degree(t, m, delta):
    n = total_num_node(t, m)
    return t + delta * n


def expect_cum_degree(x, t, m, delta):
    tx = introduce_time(x, m)
    totaldx = total_adj_degree(tx, m, delta)
    totaldt = total_adj_degree(t, m, delta)
    nnx = x + 1
    nt = total_num_node(t, m)

    b = nt - nnx
    c = totaldx / (2*m)
    d = 1 + delta/(2*m)
    prod1 = np.exp(gammaln(c/d)+gammaln(b+c/d+1/d)-gammaln(c/d+1/d)-gammaln(b+c/d))
    prod2 = totaldt / (tx + 2*m*(nt - nnx) + delta*nt)
    ksi = prod1 * prod2
    expcumd = totaldx * ksi
    return expcumd


def choose_vertex(t, eta, m, delta):
    nt = total_num_node(t, m)
    xvec = np.arange(nt)
    expcumd = expect_cum_degree(xvec, t, m, delta)
    totald = total_adj_degree(t, m, delta)
    idx = np.where((expcumd/totald) >= eta)[0]
    if not len(idx):
        raise Exception('Vertex with exp cum degree lager than {} not found!'.format(eta))
    return idx[0]


def generate_pa_edgelist(num_vertex, m, delta):
    rsd = int(time.time() * 100 % 10000)
    random.seed(rsd)
    num_tstep = introduce_time(num_vertex, m)
    tvec = range(num_tstep)
    eta_vec = [random.random() for i in range(num_tstep)]
    idx = [choose_vertex(t, eta, m, delta) for t, eta in zip(tvec, eta_vec)]
    idx = np.reshape(idx, (int(num_tstep/2), 2))
    return idx


def generate_pa_network(num_vertex, k, m, delta, max_niter=50):
    graph = nx.DiGraph()
    avg_degree = 0
    niter = 0
    while avg_degree < k and niter < max_niter:
        idx = generate_pa_edgelist(num_vertex, m, delta)
        graph.add_edges_from(idx)
        avg_degree = compute_average_degree(graph)
        print(avg_degree)
    return graph


def compute_average_degree(graph):
    return sum(dict(graph.degree()).values())/float(len(graph))
    

# def set_node_type(graph, node, node_type):
#     # Node type: 1 - excitatory, 0 - inhibitory
#     graph.nodes[node]['type'] = node_type
    

def set_excite_inhibit(graph, excite_fract):
    node_list = list(graph.nodes)
    random.shuffle(node_list)
    num_excite = int(len(node_list) * excite_fract)
    exc_list = node_list[:num_excite]
    inh_list = node_list[num_excite:]
    for u, v, d in graph.out_edges(exc_list, data=True):
        d['weight'] = 1
    for u, v, d in graph.out_edges(inh_list, data=True):
        d['weight'] = -1
    return graph


def draw_weighted_graph(G):
    eexc = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0]
    einh= [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= 0]

    pos = nx.spring_layout(G)  # positions for all nodes

    # nodes
    nx.draw_networkx_nodes(G, pos, node_size=20)

    # edges
    nx.draw_networkx_edges(G, pos, edgelist=eexc, edge_color='r', alpha=0.5)
    nx.draw_networkx_edges(G, pos, edgelist=einh, edge_color='b', alpha=0.5)

if __name__ == '__main__':
    m = 10
    delta = 1
    num_vertex = int(2.5*100)
    k = 10

    graph = generate_pa_network(num_vertex, k, m, delta)
    # graph1 = nx.read_gexf('../results/scale_free.gexf')
    graph2 = set_excite_inhibit(graph, excite_fract=0.5)

    draw_weighted_graph(graph2)
    plt.show()
    
    # nx.write_gexf(graph, "../results/scale_free.gexf")
    # nx.draw(graph, with_labels=False, node_size=20)
    # plt.show()

    # degree_sequence = sorted([d for n, d in graph.degree()], reverse=True)  # degree sequence
    # # print "Degree sequence", degree_sequence
    # degreeCount = collections.Counter(degree_sequence)
    # deg, cnt = zip(*degreeCount.items())

    # fig, ax = plt.subplots()
    # plt.plot(deg, cnt, '-.',  color='b')

    # plt.title("Degree Histogram")
    # plt.ylabel("Count")
    # plt.xlabel("Degree")
    # ax.set_xticks([d + 0.4 for d in deg])
    # ax.set_xticklabels(deg)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # plt.show()
