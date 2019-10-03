import networkx as nx
import matplotlib.pyplot as plt

DG = nx.DiGraph()
DG.add_nodes_from(range(3))
DG.add_edges_from([(1, 2), (1, 3)])


nx.draw_kamada_kawai(DG, with_labels=True, font_weight='bold')
plt.show()
