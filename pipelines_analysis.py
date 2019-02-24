import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import approximation
from nice_tree_decompositions import *

def get_graph_pipelines(shortform = True, savefig = True):
	data_file = open("./edges.csv")
	pipelines_data = data_file.readlines()
	## Process lines
	lines = []
	for line in pipelines_data[1:]:
		lines.append(line.split(";"))
	
	pipelines_data= lines
	G = nx.Graph()
	### Adding the nodes
	nodes = set([])
	for line in pipelines_data:
		if line[2] == '"Liquefied Natural Gas"':
			nodes.add('"LNG"')
		elif line[2] == '"Macedonia. Former Yugoslav Republic"':
			nodes.add('"Macedonia"')
		else:
			nodes.add(line[2])
		if line[3] == '"Liquefied Natural Gas"':
			nodes.add('"LNG"')
		elif line[3] == '"Macedonia. Former Yugoslav Republic"':
			nodes.add('"Macedonia"')
		else:
			nodes.add(line[3])
	for node in nodes:
		G.add_node(node)
	## Adding the edges
	for line in pipelines_data:
		if line[2] == '"Liquefied Natural Gas"':
			G.add_edge('"LNG"', line[3])
		elif line[3] == '"Liquefied Natural Gas"':
			G.add_edge(line[2], '"LNG"')
		elif line[2] == '"Macedonia. Former Yugoslav Republic"':
			G.add_edge('"Macedonia"', line[3])
		elif line[3] == '"Macedonia. Former Yugoslav Republic"':
			G.add_edge(line[2], '"Macedonia"')
		else:
			G.add_edge(line[2], line[3])
	# Drawing the graph
	if savefig:
		pos = nx.spring_layout(G, scale = 5)
		nx.draw(G, pos, with_labels = True, font_size = 6)
		plt.savefig("./gas_network.svg", format="svg")
		plt.clf()
	if shortform:
		G = G.subgraph(['"France"', '"Germany"', '"LNG"','"Netherlands"', '"Ukraine"', '"Switzerland"', '"Austria"', '"Italy"', '"Poland"'])
		if savefig:
			## Drawing the graph
			pos = nx.spring_layout(G, scale = 5)
			nx.draw(G, pos, with_labels = True, font_size = 6)
			plt.savefig("./small_fig.svg", format = "svg")
			plt.clf()
	return G

def main():
	G = get_graph_pipelines()
	upward_propagate = True
	edges = list(G.edges)
	edges = [frozenset(e) for e in edges]
	capacities_map = dict([])
	costs_map = dict([])
	for e in edges:
		capacities_map[e] = 1
		costs_map[e] = 2
	special_edge = edges[10]
	k = 2

	nice_tree_root_subgraphs_identifier = process_graph(G)

	special_edge = frozenset({'"Italy"', '"LNG"'})

	if upward_propagate: 
		prepared_tree_root, optimal_flow_map, min_cost, opt_solution = solve_flow(nice_tree_root_subgraphs_identifier, capacities_map, special_edge, k, costs_map)


if __name__ == "__main__":
	main()


