import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import approximation
from nice_tree_decompositions import *


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
pos = nx.spring_layout(G, scale = 5)
nx.draw(G, pos, with_labels = True, font_size = 6)
plt.savefig("./gas_network.svg", format="svg")
plt.clf()

G = G.subgraph(['"France"', '"Germany"', '"LNG"','"Netherlands"', '"Ukraine"', '"Switzerland"', '"Austria"', '"Italy"', '"Poland"'])

## Drawing the graph
pos = nx.spring_layout(G, scale = 5)
nx.draw(G, pos, with_labels = True, font_size = 6)
plt.savefig("./small_fig.svg", format = "svg")
plt.clf()




treewidth, treedecomposition = approximation.treewidth_min_fill_in(G) 

pos = nx.spring_layout(treedecomposition, scale = 3)
nx.draw(treedecomposition, pos, with_labels = True, font_size = 6)
plt.savefig("./raw_treedecomposition.svg", format="svg")
plt.clf()




edges = list(G.edges)
edges = [frozenset(e) for e in edges]
capacities_map = dict([])
costs_map = dict([])
for e in edges:
	capacities_map[e] = 1
	costs_map[e] = 2
special_edge = edges[10]
k = 3



nodes_dictionary, node_identifier_root = get_rooted_tree_decomposition(treedecomposition)
print("Tree size ", check_tree_size(nodes_dictionary[node_identifier_root], verbose = False))


nice_tree_root = get_rooted_nice_decomposition(nodes_dictionary[node_identifier_root] )
print("Nice tree size ", check_tree_size(nice_tree_root, verbose = False))
leaves =  find_leaves(nice_tree_root, verbose = False)
nice_tree_root_subgraphs, max_edges_in_subgraph = append_subgraphs(nice_tree_root, G, verbose = False)
nice_tree_root_subgraphs = clean_raw_nice_decomposition(nice_tree_root_subgraphs)



### REMOVE THISSS BELOW
tree = get_tree(nice_tree_root_subgraphs)
pos = nx.spring_layout(tree, scale = 200)
nx.draw(tree, pos, with_labels = True, font_size = 2)
plt.savefig("./raw_nice_tree_gas.svg", format ="svg")
plt.clf()





### Add all the extra edges
nice_tree_root_plus_edge = nice_tree_root_subgraphs
for special_edge in edges:
	nice_tree_root_plus_edge = add_edge_leaf(nice_tree_root_plus_edge, special_edge)
nice_tree_root_plus_edge = clean_raw_nice_decomposition(nice_tree_root_plus_edge)

nice_tree_root_subgraphs_identifier, identifier_list = append_identifier(nice_tree_root_plus_edge, verbose = False)



tree = get_tree(nice_tree_root_subgraphs_identifier)
pos = nx.spring_layout(tree, scale = 200)
nx.draw(tree, pos, with_labels = True, font_size = 2)
plt.savefig("./nice_tree_gas.svg", format ="svg")
plt.clf()


print("Find the edge")
leaves = find_leaves(nice_tree_root_subgraphs_identifier, verbose = False)
for special_edge in edges:
	found_edge = False
	for l in leaves:
		if special_edge in l.extra_info["subgraph edges"]:
			found_edge = True
	if not found_edge:
		print("Edge not found among the leaves!!!!!!!!!!!!!!")
	else:
		print("Edge found!!!")


print("Checing the stability of join nodes.")
check_join_nodes(nice_tree_root_subgraphs_identifier)
num_nice_nodes = check_tree_size(nice_tree_root_subgraphs_identifier)
print("Nice tree has ", num_nice_nodes, " nodes.")



prepared_tree_root = prepare_leaves_upward_propagation(nice_tree_root_subgraphs_identifier, 
			capacities_map ,  
			special_edge, 
			1,
			costs_map )


for i in range(num_nice_nodes):
	upward_visited = count_upward_visited(prepared_tree_root)
	print("Vertices left to upward propagate ", num_nice_nodes - upward_visited)
	prepared_tree_root = upward_propagate(prepared_tree_root,
					capacities_map ,  
					special_edge, 
					1,
					costs_map )



