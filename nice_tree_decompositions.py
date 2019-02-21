import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import approximation
import uuid
import copy 
import numpy as np
import IPython
from itertools import product
import pandas as pd

class TreeNode:
	def __init__(self, vertex_set):
		self.vertex_set = vertex_set
		self.children = []
		self.parent = None
		self.extra_info = dict([])
		self.extra_info["upward visited"] = False

	def set_children(self, tree_nodes):
		self.children = tree_nodes
	def set_parent(self, parent):
		self.parent = parent


def get_rooted_tree_decomposition( decomposition ):
	nodes_identifiers = list(decomposition.nodes) ## list of frozensets
	#edges = list(treedecomposition.edges) ## list of tuples having two frozenset
	nodes_dictionary = dict([])
	for node_identifier in nodes_identifiers:
		nodes_dictionary[node_identifier] = TreeNode(node_identifier)

	## Set the root to be the first node.
	node_identifier_root = nodes_identifiers[0]
	frontier_node_identifiers = [node_identifier_root]
	visited_node_identifiers = []

	while frontier_node_identifiers != []:
		focus_node_identifier = frontier_node_identifiers[0]
		if focus_node_identifier in visited_node_identifiers:
			frontier_node_identifiers.remove(focus_node_identifier)
			continue
		neighbor_identifiers = list(decomposition.neighbors(focus_node_identifier))
		focus_node = nodes_dictionary[focus_node_identifier]
		focus_parent = focus_node.parent
		focus_parent_identifier = None
		if focus_parent != None:
			focus_parent_identifier = focus_parent.vertex_set

		neighbor_list = []
		for neighbor_identifier in neighbor_identifiers:
			if neighbor_identifier != focus_parent_identifier:
				neighbor_list.append(nodes_dictionary[neighbor_identifier])
				nodes_dictionary[neighbor_identifier].set_parent(focus_node)
				frontier_node_identifiers.append(neighbor_identifier)


		focus_node.set_children(neighbor_list)
		frontier_node_identifiers = list(set(frontier_node_identifiers))

		visited_node_identifiers.append(focus_node_identifier)
		frontier_node_identifiers.remove(focus_node_identifier)
	return nodes_dictionary, node_identifier_root





def count_upward_visited(root_treenode, verbose = False):
	counter= 0 
	frontier = [root_treenode]
	while frontier != []:
		focus_node = frontier[0]
		if verbose:
			print(focus_node.vertex_set)
		if focus_node.extra_info["upward visited"]:
			counter += 1
		frontier = frontier[1:]
		frontier += focus_node.children
	return counter



def check_tree_size(root_treenode, verbose = False):
	counter= 0 
	frontier = [root_treenode]
	while frontier != []:
		focus_node = frontier[0]
		if verbose:
			print(focus_node.vertex_set)
		counter += 1
		frontier = frontier[1:]
		frontier += focus_node.children
	return counter

def find_leaves(root_treenode, verbose = False):
	leaves = []
	frontier = [root_treenode]
	while frontier != []:
		focus_node = frontier[0]
		if verbose:
			print(focus_node.vertex_set)
		if len(focus_node.children) == 0:
			leaves.append(focus_node)
		frontier = frontier[1:]
		frontier += focus_node.children
	return leaves

def append_subgraphs(root_treenode, G, verbose = True):
	tree_with_subgraphs = copy.deepcopy(root_treenode)
	max_edges_in_subgraph = -float("inf")
	frontier = [tree_with_subgraphs]
	while frontier != []:
		focus_node = frontier[0]
		focus_node.extra_info["subgraph"] = G.subgraph(list(focus_node.vertex_set))
		subgraph_edges = list(focus_node.extra_info["subgraph"].edges)
		subgraph_edges.sort()
		subgraph_edges = [frozenset(e) for e in subgraph_edges]
		focus_node.extra_info["subgraph edges"] = subgraph_edges 		

		max_edges_in_subgraph = max(max_edges_in_subgraph, len(focus_node.extra_info["subgraph"].edges))
		if verbose:
			print(focus_node.vertex_set)
			print("Max edges in subraph so far ", max_edges_in_subgraph)
		frontier = frontier[1:]
		frontier += focus_node.children
	return tree_with_subgraphs, max_edges_in_subgraph

def append_identifier(root_treenode, verbose = True):
	tree_with_subgraphs = copy.deepcopy(root_treenode)
	identifier_list = []
	frontier = [tree_with_subgraphs]
	while frontier != []:
		focus_node = frontier[0]
		focus_node.extra_info["identifier"] = uuid.uuid4() 
		identifier_list.append(focus_node.extra_info["identifier"])
		if verbose:
			print(focus_node.vertex_set)
		frontier = frontier[1:]
		frontier += focus_node.children
	return tree_with_subgraphs, identifier_list

def check_edge_belonging(root_treenode, target_edge, verbose = True):
	frontier = [root_treenode]
	while frontier != []:
		focus_node = frontier[0]
		subgraph_edges = focus_node.extra_info["subgraph edges"]
		if target_edge in subgraph_edges:
			print("Target edge found!!!!")

		if verbose:
			print(focus_node.vertex_set)
		frontier = frontier[1:]
		frontier += focus_node.children

def check_join_nodes(root_treenode, verbose = False):
	frontier = [root_treenode]
	while frontier != []:
		focus_node = frontier[0]
		if len(focus_node.children) == 2:
			### Perform checks
			if len(focus_node.extra_info["subgraph edges"]) != len(focus_node.children[0].extra_info["subgraph edges"]):
				raise ValueError("This join node is wrong. One of its children has different number of edges.")
			if len(focus_node.extra_info["subgraph edges"]) != len(focus_node.children[1].extra_info["subgraph edges"]):
				raise ValueError("This join node is wrong. One of its children has different number of edges.") 

			for i in range(len(focus_node.extra_info["subgraph edges"])):
				if focus_node.extra_info["subgraph edges"][i] != focus_node.children[0].extra_info["subgraph edges"][i]:
					raise ValueError("This join node is wrong! Edge sets do not match")
				if focus_node.extra_info["subgraph edges"][i] != focus_node.children[1].extra_info["subgraph edges"][i]:
					raise ValueError("This join node is wrong! Edge sets do not match")
			focus_vertex_set = list(focus_node.vertex_set)
			child_1_vertex_set = list(focus_node.children[0].vertex_set)
			child_2_vertex_set = list(focus_node.children[1].vertex_set)
			if len(focus_vertex_set ) != len(child_1_vertex_set):
				raise ValueError("This join node is wrong. One of its children has different number of vertices.")
			if len(focus_vertex_set ) != len(child_2_vertex_set):
				raise ValueError("This join node is wrong. One of its children has different number of vertices.")
			for j in range(len(focus_node.vertex_set)):  
				if focus_vertex_set[j] != child_1_vertex_set[j]:
					raise ValueError("This join node is wrong!. Vertex sets do not match.")
				if focus_vertex_set[j] != child_2_vertex_set[j]:
					raise ValueError("This join node is wrong!. Vertex sets do not match.")
		frontier = frontier[1:]
		frontier += focus_node.children

 


def get_rooted_nice_decomposition( root_treenode ):
	nice_tree = copy.deepcopy(root_treenode)
	frontier = [nice_tree]
	while frontier != []:
		focus_node = frontier[0]
		frontier = frontier[1:]
		focus_children = focus_node.children

		### If focus_node is a leaf node.
		if len(focus_children) == 0:
			continue

		focus_child = focus_children[0]
		other_children = focus_children[1:]

		if len(other_children) > 0:
			focus_node_join_1 = TreeNode(focus_node.vertex_set)
			focus_node_join_1.parent = focus_node
			focus_node_join_1.children = [focus_child]
			focus_child.parent = focus_node_join_1

			focus_node_join_2 = TreeNode(focus_node.vertex_set)
			focus_node_join_2.parent = focus_node
			focus_node_join_2.children = other_children
			for child in other_children:
				child.parent = focus_node_join_2

			focus_node.children = [focus_node_join_1, focus_node_join_2]
			frontier.append(focus_node_join_2)

			focus_node = focus_node_join_1
		### Find the necessary set operations to turn the identifier of focus node into the 
		### identifier of the focus child
		
		focus_identifier_set = set(focus_node.vertex_set)
		focus_child_identifier_set = set(focus_child.vertex_set)

		elements_to_remove = focus_identifier_set - focus_child_identifier_set
		elements_to_add = focus_child_identifier_set - focus_identifier_set

		current_identifier_set = copy.copy(focus_identifier_set)
		identifier_sets_in_chain = []

		for e in elements_to_remove:
			current_identifier_set.remove(e)
			identifier_sets_in_chain.append(frozenset(current_identifier_set))
		for e in elements_to_add:
			current_identifier_set.add(e)
			identifier_sets_in_chain.append(frozenset(current_identifier_set))

		### Check that the last set in the chain equals the set in focus_cchild
		if set(identifier_sets_in_chain[-1]) != set(focus_child_identifier_set):
			raise ValueError("The chain of additions and subtractions was not right.")

		identifier_sets_in_chain = identifier_sets_in_chain[:-1]

		### In case we have to add a chain, create the nodes and redo the connections.
		if len(identifier_sets_in_chain) > 0:
			chain_nodes = []

			### Create the chain nodes
			for identifier_set in identifier_sets_in_chain:
				chain_nodes.append(TreeNode(identifier_set))

			### Create the connections
			focus_node.children = [chain_nodes[0]]

			for i in range(len(chain_nodes)):
				if i ==0:
					chain_nodes[i].parent = focus_node
				else:
					chain_nodes[i].parent = chain_nodes[i-1]

				if i == len(chain_nodes)-1:
					chain_nodes[i].children = [focus_child]
				else:
					chain_nodes[i].children = [chain_nodes[i+1]]

			focus_child.parent = chain_nodes[-1]

		frontier.append(focus_child)

	return nice_tree

## compresses the graph so it doesn't have node parent structures where the parent has a single 
## child and the child has the same vertex set as the parent.
def clean_raw_nice_decomposition( root_treenode ):
	nice_tree = copy.deepcopy(root_treenode)
	frontier = [nice_tree]
	while frontier != []:
		focus_node = frontier[0]
		#print("Child alsdkmfalsdmalksdmfalksdmflaksdmlaksmdflkm ",focus_node.vertex_set, [c.vertex_set for c in focus_node.children] )
		
		if len(focus_node.children) == 1:
			#print("Child of size 1")
			focus_child = focus_node.children[0]
			if focus_child.vertex_set == focus_node.vertex_set:
				#print("Cleaning child parent structure")
				focus_node.children = focus_child.children
				focus_child.parent = None
				for focs_granchild in focus_child.children:
					focus_granchild.parent = focus_node
				focus_child.children = None

		frontier = frontier[1:]
		frontier += focus_node.children
	return nice_tree


def add_edge_leaf(root_treenode, target_edge):
	root_edge_added = copy.deepcopy(root_treenode)
	leaves = find_leaves(root_edge_added, verbose = False)
	for l in leaves:
		if target_edge in l.extra_info["subgraph edges"]:
			print("Edge was already there!!!")
			return root_edge_added

	frontier = [root_edge_added]
	while frontier != []:
		focus_node = frontier[0]
		subgraph_edges = focus_node.extra_info["subgraph edges"]
		if target_edge in subgraph_edges:
			edge_found = True
			break
		frontier = frontier[1:]
		frontier += focus_node.children
	if not edge_found:
		raise ValueError("Edge is not found!!!!!")

	focus_subgraph = focus_node.extra_info["subgraph"]
	focus_subgraph_edges = focus_node.extra_info["subgraph edges"] 

	if len(focus_node.children) == 1:
		focus_child = focus_node.children[0]
		if focus_child.vertex_set == focus_node.vertex_set:
			focus_node_join_1 = focus_child
		else:
			focus_node_join_1 = TreeNode(focus_node.vertex_set)
			focus_node_join_1.children = [focus_child]
			focus_node_join_1.parent = focus_node
			focus_child.parent = focus_node_join_1

		focus_node_join_2 = TreeNode(focus_node.vertex_set)
		focus_node_join_2.parent = focus_node
		focus_node.children = [focus_node_join_1, focus_node_join_2]
		focus_node = focus_node_join_2

	if len(focus_node.children) == 2:
		focus_children = focus_node.children
		focus_node_join_1 = TreeNode(focus_node.vertex_set)
		focus_node_join_1.parent = focus_node
		focus_node_join_1.children = focus_children
		focus_children[0].parent = focus_node_join_1
		focus_children[1].parent = focus_node_join_1
		focus_node_join_2 = TreeNode(focus_node.vertex_set)
		focus_node_join_2.parent = focus_node
		focus_node.children = [focus_node_join_1, focus_node_join_2]
		focus_node = focus_node_join_2

	focus_node_join_1.extra_info["subgraph"] = copy.deepcopy(focus_subgraph)
	focus_node_join_1.extra_info["subgraph edges"] = copy.copy(focus_subgraph_edges)
	focus_node_join_2.extra_info["subgraph"] = copy.deepcopy(focus_subgraph)
	focus_node_join_2.extra_info["subgraph edges"] = copy.copy(focus_subgraph_edges)

	edge_leaf = TreeNode(frozenset(target_edge))
	edge_leaf.parent = focus_node
	edge_leaf.extra_info["subgraph"] = focus_subgraph.subgraph(target_edge)

	edge_leaf_edges = list(edge_leaf.extra_info["subgraph"].edges)
	edge_leaf_edges = [frozenset(e) for e in edge_leaf_edges]

	edge_leaf.extra_info["subgraph edges"] = edge_leaf_edges

	focus_node.children = [edge_leaf]

	focus_identifier_set = set(focus_node.vertex_set)
	edge_leaf_identifier_set = set(target_edge)

	elements_to_remove = focus_identifier_set - edge_leaf_identifier_set
	elements_to_add = edge_leaf_identifier_set - focus_identifier_set

	current_identifier_set = copy.copy(focus_identifier_set)
	identifier_sets_in_chain = []
	for e in elements_to_remove:
		current_identifier_set.remove(e)
		identifier_sets_in_chain.append(frozenset(current_identifier_set))
	for e in elements_to_add:
		current_identifier_set.add(e)
		identifier_sets_in_chain.append(frozenset(current_identifier_set))

	#print("current indentifier set ", current_identifier_set)
	#print("edge leaf identifier set ", edge_leaf_identifier_set)
	#print(identifier_sets_in_chain)
	### check that the last set in the chain equals the edge
	if len(focus_node.vertex_set) > 2:
		if set(identifier_sets_in_chain[-1]) != set(edge_leaf_identifier_set):
			raise ValueError("The chain of additions and subtractions was not right")

	identifier_sets_in_chain = identifier_sets_in_chain[:-1]



	### Create the chain nodes and redo the connections
	if len(identifier_sets_in_chain) > 0:
		chain_nodes = []

		### Create the chain nodes
		for identifier_set in identifier_sets_in_chain:
			chain_nodes.append(TreeNode(identifier_set))
			chain_nodes[-1].extra_info["subgraph"] = focus_subgraph.subgraph(identifier_set)
			chain_node_edges = list(chain_nodes[-1].extra_info["subgraph"].edges)
			chain_node_edges = [frozenset(e) for e in chain_node_edges]
			chain_nodes[-1].extra_info["subgraph edges"]  = chain_node_edges

		focus_node.children = [chain_nodes[0]]

		for i in range(len(chain_nodes)):
			if i == 0:
				chain_nodes[i].parent = focus_node
			else:
				chain_nodes[i].parent = chain_nodes[i-1]

			if i == len(chain_nodes)-1:
				chain_nodes[i].children = [edge_leaf]
			else:
				chain_nodes[i].children = [chain_nodes[i+1]]
		edge_leaf.parent = chain_nodes[-1]

	return root_edge_added



def add_capacities_and_costs(root_treenode, capacities_map, costs_map, verbose = False):
	tree_with_costs_capacities = copy.deepcopy(root_treenode)
	frontier = [tree_with_costs_capacities]
	while frontier != []:
		focus_node = frontier[0]

		if verbose:
			print(focus_node.vertex_set)
		frontier = frontier[1:]
		frontier += focus_node.children
	return tree_with_subgraphs


def get_flow_into_vertex(frozenset_edge, target_vertex, flow_value):
	edge_list = list(frozenset_edge)
	#edge_list[0] ---> edge_list[1]
	# The flow is assumed to go from the first vertex to the second one.
	if edge_list[0] == target_vertex:
		return - flow_value
	elif edge_list[1] == target_vertex:
		return flow_value
	else:
		raise ValueError("Flow into vertex function's frozenset edge doesn't have the target vertex.")

def get_tree(root_treenode, verbose = False):
	frontier = [root_treenode]
	counter = 1
	root_treenode.extra_info["tree counter"] = counter
	while frontier != []:
		focus_node = frontier[0]
		focus_node.extra_info["tree counter"] = counter
		counter += 1
		if verbose:
			print(focus_node.vertex_set)
		frontier = frontier[1:]
		frontier += focus_node.children


	frontier = [root_treenode]
	tree  = nx.Graph()
	counter = 1
	while frontier != []:
		focus_node = frontier[0]
		for child in focus_node.children:
			tree.add_edge((tuple(focus_node.vertex_set), focus_node.extra_info["tree counter"]), (tuple(child.vertex_set), child.extra_info["tree counter"]))

		if verbose:
			print(focus_node.vertex_set)
		frontier = frontier[1:]
		frontier += focus_node.children
	return tree


def format_edge_name(edge):
	edge_as_tuple = tuple(edge)
	edge_name = edge_as_tuple[0]+ "-" + edge_as_tuple[1]
	return edge_name

def prepare_leaves_upward_propagation(root_treenode, capacities_map, special_edge,  k, costs_map):
	root_nice_treenode_subgraph = copy.deepcopy(root_treenode)
	special_edge = frozenset(special_edge)
	root_nice_treenode_subgraph_identifier, identifier_list = append_identifier(root_nice_treenode_subgraph, verbose = False)
	
	leaves= find_leaves(root_nice_treenode_subgraph_identifier, verbose = False)
	## Initialize flows at leaves
	for leaf in leaves:
		subgraph_edges = leaf.extra_info["subgraph edges"]
		edges_to_indices = dict([(subgraph_edges[i], i) for i in range(len(subgraph_edges))])
		subgraph_capacities = [capacities_map[subgraph_edges[i]] for i in range(len(subgraph_edges))]
		subgraph_capacities_list = [list(range(-c-1, c+1)) for c in subgraph_capacities]
		leaf_vertices = list(leaf.vertex_set)
		vertices_to_indices = dict([(leaf_vertices[i], i) for i in range(len(leaf_vertices))])
		leaf.extra_info["edges to indices"] = edges_to_indices 
		leaf.extra_info["vertices to indices"] = vertices_to_indices

		#df_columns = leaf_vertices + [format_edge_name(e) for e in subgraph_edges] + ["cost"] 
		#partial_solutions_df = pd.DataFrame( columns = df_columns)

		partial_solutions = []
		if special_edge in subgraph_edges:
			print("Special edge is in this leaf!!")
			special_edge_index = edges_to_indices[special_edge]

			subgraph_capacities_list[special_edge_index] = [k,]
		local_flows = product(*subgraph_capacities_list)
		leaf.extra_info["local flows"] = local_flows
		leaf.extra_info["subgraph capacities list"] = subgraph_capacities_list

		#partial_solutions = []
		for frontier_flow in local_flows:
			residual_flow = np.zeros(len(leaf_vertices))
			x_cost = 0
			for edge in leaf.extra_info["subgraph edges"]:
				edge_index = edges_to_indices[edge]
				edge_flow = frontier_flow[edge_index]
				if edge_flow != 0:
					x_cost += costs_map[edge]
				
			for i in range(len(leaf_vertices)):
				v = leaf_vertices[i]
				neighbors = list(leaf.extra_info["subgraph"].neighbors(v))
				for neighbor in neighbors:
					edge = frozenset([neighbor, v])
					edge_index = edges_to_indices[edge]
					edge_flow = frontier_flow[edge_index]
					if edge_flow != 0:
						residual_flow[i] +=  get_flow_into_vertex(edge, v, edge_flow)
			
			partial_solutions.append((residual_flow, frontier_flow, x_cost))
			#dataframe_row = list(residual_flow) + list(frontier_flow) + [x_cost]
			#row_df = pd.DataFrame(dataframe_row, columns = df_columns)
			#partial_solutions_df.append(row_df, ignore_index = True)


		leaf.extra_info["partial solutions"] = partial_solutions
		leaf.extra_info["upward visited"] = True

	return root_nice_treenode_subgraph_identifier

def upward_propagate(root_treenode, capacities_map, special_edge,  k, costs_map, verbose = False):
	frontier = [root_treenode]
	already_finished = True
	while frontier != []:
		focus_node = frontier[0]
		subgraph_edges = focus_node.extra_info["subgraph edges"]
		if special_edge in subgraph_edges and verbose:
			print("Target edge found!!!!")

		visited_focus_node = focus_node.extra_info["upward visited"]
		visited_children = [c.extra_info["upward visited"] for c in focus_node.children]

		if not visited_focus_node and np.sum(visited_children) == len(visited_children):
			already_finished = False
			break
		frontier = frontier[1:]
		frontier += focus_node.children
	
	if already_finished:
		print("Already finished")
		return root_treenode


	subgraph_edges = focus_node.extra_info["subgraph edges"]
	edges_to_indices = dict([(subgraph_edges[i], i) for i in range(len(subgraph_edges))])
	subgraph_capacities = [capacities_map[subgraph_edges[i]] for i in range(len(subgraph_edges))]
	subgraph_capacities_list = [list(range(-c-1, c+1)) for c in subgraph_capacities]

	focus_node.extra_info["edges to indices"] = edges_to_indices
	focus_node_vertices = list(focus_node.vertex_set)
	vertices_to_indices = dict([(focus_node_vertices[i], i) for i in range(len(focus_node_vertices))])
	focus_node.extra_info["vertices to indices"] = vertices_to_indices


	#### FIXING THE SUBGRAPH CAPACITIES LIST TO ACCOUNT FOR THE SPECIAL EDGE
	subgraph_vertices = list(focus_node.vertex_set)
	vertices_to_indices = dict([(subgraph_vertices[i], i) for i in range(len(subgraph_vertices))])

	if special_edge in subgraph_edges:
		print("Special edge is in this node!!")
		special_edge_index = edges_to_indices[special_edge]
		subgraph_capacities_list[special_edge_index] = [k,]


	focus_vertices = list(focus_node.vertex_set)
	vertices_to_indices = dict([(focus_vertices[i], i) for i in range(len(focus_vertices))])
	#df_columns = subgraph_vertices + [format_edge_name(e) for e in subgraph_edges] + ["cost"] 
	#partial_solutions_df = pd.DataFrame(columns = df_columns )
	partial_solutions = []

	if len(focus_node.children) == 2 and len(focus_node.vertex_set) == len(focus_node.children[0].vertex_set) and len(focus_node.vertex_set) == len(focus_node.children[0].vertex_set):
		print("Found a join node!!!")

		print("Children partial solution tables sizes ", len(focus_node.children[0].extra_info["partial solutions"]), len(focus_node.children[1].extra_info["partial solutions"]))
		## CREATE DICTIONARY FROM THE PARTIAL SOLUTIONS LISTS
		partial_solutions_map_1 = dict([])
		for partial_solution in focus_node.children[0].extra_info["partial solutions"]:
			(child_residual_flow_1, child_frontier_flow_1, child_x_cost_1) = partial_solution 
			partial_solution_key = tuple( child_frontier_flow_1 )
			if partial_solution_key not in partial_solutions_map_1:
				partial_solutions_map_1[partial_solution_key] = [( child_residual_flow_1, child_x_cost_1)]
			else:
				partial_solutions_map_1[partial_solution_key].append(( child_residual_flow_1, child_x_cost_1))


		for partial_solution_2 in focus_node.children[1].extra_info["partial solutions"]:
			(child_residual_flow_2, child_frontier_flow_2, child_x_cost_2) = partial_solution_2 
			partial_solution_key = tuple( child_frontier_flow_2 )
			
			if partial_solution_key not in partial_solutions_map_1:
				continue

			child_1_matching_solutions = partial_solutions_map_1[partial_solution_key]

			for (child_residual_flow_1, child_x_cost_1_retrieved) in child_1_matching_solutions:

				joined_cost = child_x_cost_1_retrieved + child_x_cost_2
				frontier_flow_cost_component = 0

				## Subtract the double counting of costs from the frontier flow
				child_edges  = focus_node.children[0].extra_info["subgraph edges"]
				for j in range(len(child_edges)):
					child_edge = child_edges[j]
					if child_frontier_flow_2[j] > 0:
						joined_cost -= costs_map[child_edge]

				frontier_flow = child_frontier_flow_2 
				x_cost = joined_cost
				### Compute the new resulting residual flow.
				residual_flow = list( np.array( child_residual_flow_1) + np.array(child_residual_flow_2) )


				for i in range(len(focus_vertices)):
					v = focus_vertices[i]
					neighbors = list(focus_node.extra_info["subgraph"].neighbors(v))
					for neighbor in neighbors:
						edge = frozenset([neighbor, v])
						edge_index = edges_to_indices[edge]
						edge_flow = frontier_flow[edge_index]
						if edge_flow != 0:
							residual_flow[i] -=  get_flow_into_vertex(edge, v, edge_flow)

				partial_solutions.append((residual_flow, frontier_flow, x_cost))




	elif len(focus_node.children) == 1 and len(focus_node.vertex_set) == len(focus_node.children[0].vertex_set) -1:
		print("Found a forget node!!")
		print("Children partial solution tables sizes ", len(focus_node.children[0].extra_info["partial solutions"]) )
		forgotten_vertex = list(focus_node.children[0].vertex_set - focus_node.vertex_set)[0]
		neighbors = list(focus_node.children[0].extra_info["subgraph"].neighbors(forgotten_vertex))
		child_edges_to_indices = focus_node.children[0].extra_info["edges to indices"]
		child_edges = focus_node.children[0].extra_info["subgraph edges"]
		child_vertices_to_indices = focus_node.children[0].extra_info["vertices to indices"]

	
		projected_solutions_map = dict([])
		for partial_solution in focus_node.children[0].extra_info["partial solutions"]:
			(child_residual_flow, child_frontier_flow, child_x_cost) = partial_solution
			
			## Check if residual flow at the forgotten vertex is zero
			if child_residual_flow[child_vertices_to_indices[forgotten_vertex]] != 0:
				continue



			projected_residual_flow = [child_residual_flow[child_vertices_to_indices[v]] for v in focus_vertices]
			projected_partial_flow = [child_frontier_flow[child_edges_to_indices[e]] for e in subgraph_edges]
			neighbors_cost = 0

			# for n in neighbors:
			# 	neighbor_edge = frozenset([n, forgotten_vertex])
			# 	neighbor_edge_index = child_edges_to_indices[neighbor_edge]
			# 	neighbor_edge_flow = child_residual_flow[neighbor_edge_index]
			# 	if neighbor_edge_flow != 0:
			# 		neighbors_cost += costs_map[neighbor_edge]



			projected_cost = child_x_cost + neighbors_cost ## I don't think we need to modify the cost value.

			projected_solution_key = tuple( projected_residual_flow + projected_partial_flow )     ## Hash the residual flow profile and the 

			if projected_solution_key in projected_solutions_map:
				projected_solutions_map[projected_solution_key] = min( projected_solutions_map[projected_solution_key],  projected_cost) 
			else:
				projected_solutions_map[projected_solution_key] = projected_cost

		for projected_solution_key in projected_solutions_map:
			residual_flow =list( projected_solution_key[ :len(subgraph_vertices) ] ) 
			frontier_flow = list(projected_solution_key[len(subgraph_vertices): len(subgraph_vertices) + len(subgraph_edges) ])
			x_cost = projected_solutions_map[projected_solution_key]

			partial_solutions.append((residual_flow, frontier_flow, x_cost))





	elif len(focus_node.children) == 1 and len(focus_node.vertex_set) == len(focus_node.children[0].vertex_set) +1:
		print("Found an introduce node!!")
		print("Children partial solution tables sizes ", len(focus_node.children[0].extra_info["partial solutions"]) )

		added_vertex = list(focus_node.vertex_set - focus_node.children[0].vertex_set)[0]
		added_vertex_index = focus_node.extra_info["vertices to indices"][added_vertex]
		neighbors = list(focus_node.extra_info["subgraph"].neighbors(added_vertex))
		child_edges_to_indices = focus_node.children[0].extra_info["edges to indices"]
		child_vertices_to_indices = focus_node.children[0].extra_info["vertices to indices"]

		edges_neighbors = [frozenset([n, added_vertex]) for n in neighbors]
		edges_neighbors_indices = [edges_to_indices[e] for e in edges_neighbors]
		neighbors_capacities_list = [subgraph_capacities_list[neighbor_edge_index] for neighbor_edge_index in edges_neighbors_indices]


		for neighbor_flow_profile in product(*neighbors_capacities_list):
			
			for partial_solution in focus_node.children[0].extra_info["partial solutions"]: 
				#partial_solution = numpy.array(partial_solution)
				#if len(partial_solution) !=  len(child_vertices_to_indices) + len(child_edges_to_indices) + 1:
				#	raise ValueError("The child partial solution in this introduce node does not match the expected size")
				(child_residual_flow, child_frontier_flow, child_x_cost) = partial_solution
				#child_residual_flow = partial_solution[:len(child_vertices_to_indices)]
				#child_frontier_flow = partial_solution[len(child_vertices_to_indices): len(child_vertices_to_indices) + len(child_edges_to_indices)]
				#child_x_cost = partial_solution[-1]
				(child_residual_flow, child_frontier_flow, child_x_cost) = partial_solution
				residual_flow = np.zeros(len(focus_node.vertex_set))
				frontier_flow = np.zeros(len(focus_node.extra_info["subgraph edges"]))
				x_cost = partial_solution[2]
				
				for edge in focus_node.children[0].extra_info["subgraph edges"]:
					child_edge_index = child_edges_to_indices[edge]
					edge_index = edges_to_indices[edge]

					frontier_flow[edge_index] = child_frontier_flow[child_edge_index]

				for vertex in focus_node.children[0].extra_info["vertices to indices"]:
					child_vertex_index = child_vertices_to_indices[vertex]
					vertex_index = vertices_to_indices[vertex]
					residual_flow[vertex_index] = child_residual_flow[child_vertex_index]


				for j in range(len(neighbors)):
					n = neighbors[j]
					edge = frozenset([n, added_vertex])
					edge_index = edges_to_indices[edge]
					edge_flow = neighbor_flow_profile[j]

					frontier_flow[edge_index] = edge_flow
					neighbor_vertex_index = vertices_to_indices[n]
					residual_flow[added_vertex_index] += get_flow_into_vertex(edge, added_vertex, edge_flow)
					residual_flow[neighbor_vertex_index] += get_flow_into_vertex(edge, n, edge_flow)
					if edge_flow != 0:
						x_cost += costs_map[edge]

				#dataframe_row = list(residual_flow) + list(frontier_flow) + [x_cost]
				#row_df = pd.DataFrame(dataframe_row, columns = df_columns )
				#partial_solutions_df.append(row_df, ignore_index = True)
				partial_solutions.append((residual_flow, frontier_flow, x_cost))

	else:
		raise ValueError("This parent - child structure is broken!!!")
	#focus_node.extra_info["partial solutions"] = partial_solutions_df
	focus_node.extra_info["upward visited"] = True
	focus_node.extra_info["partial solutions"] = partial_solutions
	print("Size of the partial solutions ", len(partial_solutions))
	return root_treenode



