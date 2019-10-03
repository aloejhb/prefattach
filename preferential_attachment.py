import time
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.special import gamma
from operator import mul
from functools import reduce


def introduce_time(x, m):
    return x * 2*m


def total_num_node(t, m):
    return np.floor(t/(2*m)) + 1


def total_adj_degree(t, m, delta):
    n = total_num_node(t, m)
    return t + delta * n


def get_ksi(c, d, istep, m, delta):
    plist = [1+1/(c+i*d) for i in range(istep)]
    return reduce(mul, plist, 1)
    

def expect_cum_degree(x, t, m, delta):
    tx = introduce_time(x, m)
    totaldx = total_adj_degree(tx, m, delta)
    totaldt = total_adj_degree(t, m, delta)
    nnx = total_num_node(tx, m)
    nt = total_num_node(t, m)
    b = nt - nnx
    if b == 0:
        ksi = (t + delta*nnx) / (tx + delta*nnx)
    else:
        c = totaldx / (2*m)
        d = 1 + delta / (2*m)
        prod1 = get_ksi_second_prod(istep, c, m, delta)
        prod2 = totaldt / (tx+2*m*b+delta*nt)
        ksi = prod1 * prod2
    expcumd =  totaldx * ksi
    if np.isnan(expcumd).any():
        import pdb; pdb.set_trace()
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


def generate_pa_network(num_vertex, m, delta):
    np.random.seed(int(time.time()))
    num_tstep = introduce_time(num_vertex, m)
    tvec = range(num_tstep)
    eta_vec = np.random.rand(num_tstep)
    idx = [choose_vertex(t, eta, m, delta) for t, eta in zip(tvec, eta_vec)]
    idx = np.reshape(idx, (int(num_tstep/2), 2))
    graph = nx.from_edgelist(idx, create_using=nx.MultiDiGraph)
    return graph
    


if __name__ == '__main__':
    m = 10
    delta = 1
    num_vertex = 200

    tx = 2000
    t = 4000
    nnx = total_num_node(tx, m)
    nt = total_num_node(t, m)
    
    c = (tx + delta*nnx)/(2*m)
    d = 1+delta/(2*m)
    istep = int(nt-nnx)
    ans = get_ksi(c, d, istep, m, delta)
    # graph = generate_pa_network(num_vertex, m, delta)
    
    # nx.draw_kamada_kawai(graph, with_labels=True, font_weight='bold')
    # plt.show()

