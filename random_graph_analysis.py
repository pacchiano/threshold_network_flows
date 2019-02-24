import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import approximation
from nice_tree_decompositions import *

def get_geometric_graph(n, radius):
	G = nx.random_geometric_graph(n, radius)
	return G
def get_duplication_graph(n, p, special_edge = None):
	G = nx.duplication_divergence_graph(n, p)
	if special_edge != None:
		G.add_edge(special_edge[0], special_edge[1])
	return G

def main():
	G = get_geometric_graph(n=20, radius=.15)
	print("Graph has ", len(G.edges), " edges.")
	print("Graph has ", len(G.nodes), " vertices.")
	upward_propagate = True
	special_edge_index = 10
	edges = list(G.edges)
	edges = [frozenset(e) for e in edges]
	capacities_map = dict([])
	costs_map = dict([])
	for e in edges:
		capacities_map[e] = 1
		costs_map[e] = 2
	k = 1
	#IPython.embed()
	nice_tree_root_subgraphs_identifier = process_graph(G)
	special_edge = edges[min(special_edge_index, len(G.edges)-1)]

	if upward_propagate: 
		prepared_tree_root, optimal_flow_map, min_cost, opt_solution = solve_flow(nice_tree_root_subgraphs_identifier, capacities_map, special_edge, k, costs_map)

	print("Min cost ", min_cost)


if __name__ == "__main__":
	main()