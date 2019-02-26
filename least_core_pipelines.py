import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from scipy import optimize 
from networkx.algorithms import approximation
from nice_tree_decompositions import *
from pipelines_analysis import *
from random_graph_analysis import *




def get_epsilon_core( k = 1, use_tree = True, 
	special_edge = None, G = None, prunning_cost = float("inf"), max_num_constraints = float("inf")):

	if not use_tree and k > 1:
		raise ValueError("Not using tree and k > 1")
	epsilon_infeasible = 0
	epsilon_feasible = 5
	tolerance = .05
	edges = list(G.edges)

	vertices = list(G.nodes)
	edges = [frozenset(e) for e in edges]
	edges_to_indices = dict([])
	capacities_map = dict([])
	for i in range(len(edges)):
		e = edges[i]
		capacities_map[e] = 1
		edges_to_indices[e] = i 

	#IPython.embed()
	#print(edges_to_indices)
	special_edge_index = edges_to_indices[special_edge]
	coalitions = []
	imputation = np.ones(len(edges))*1.0/(len(edges)-1)
	imputation[special_edge_index] = 0

	if use_tree:
		nice_tree_root_subgraphs_identifier = process_graph(G)

	while epsilon_feasible - epsilon_infeasible > tolerance:
		costs_map = dict([])
		for e in edges:
			index = edges_to_indices[e]
			costs_map[e] = imputation[index]
			print("Edge ", e, " imputation value ", imputation[index])
		if use_tree:
			prepared_tree_root, optimal_flow_map, min_cost, opt_solution = solve_flow(nice_tree_root_subgraphs_identifier, 
												capacities_map, special_edge, k, costs_map)
		else:
			cost_vector = np.zeros(len(edges))
			for edge in edges:
				index = edges_to_indices[edge]
				cost_vector[index]= imputation[index]

			flow_constraints = np.zeros((len(vertices), len(edges)))
			flow_equals = np.zeros(len(vertices))
			for i in range(len(vertices)):
				v = vertices[i]
				neighbors = list(G.neighbors(v))
				for n in neighbors:
					edge = frozenset([v,n])
					index = edges_to_indices[edge]
					flow_constraints[i,index] = get_flow_into_vertex(edge, v, 1)

			capacity_bounds = np.ones(2*len(edges))
			capacity_constraints = np.vstack((np.eye(len(edges)) , -np.eye(len(edges))))
			special_edge_constraint = np.zeros(len(edges))
			special_edge_constraint[special_edge_index] = 1
			flow_constraints = np.vstack((flow_constraints, special_edge_constraint))
			flow_equals = np.zeros(len(vertices) +1)
			flow_equals[-1] = k
			bounds = list(zip(-np.ones(len(edges)),np.ones(len(edges))))

			### we add absolute value variables
			abs_constraints_top = np.vstack((np.eye(len(edges)) , -np.eye(len(edges))))
			abs_constraints_bottom =  np.vstack((-np.eye(len(edges)) , -np.eye(len(edges))))
			abs_constraints = np.hstack((abs_constraints_top, abs_constraints_bottom))
			abs_constraints_right = np.zeros(2*len(edges))
			cost_vector = np.hstack((np.zeros(len(edges)), cost_vector))
			bounds += list(zip(len(edges)*[0], len(edges)*[None]))

			flow_constraints = np.hstack((flow_constraints, np.zeros((len(vertices) +1, len(edges)))))

			print("Starting the flow optimization!")
			opt_flow = optimize.linprog(cost_vector, A_eq= flow_constraints,
			 b_eq = flow_equals, bounds = bounds, A_ub = abs_constraints, b_ub = abs_constraints_right )#A_ub =capacity_constraints, b_ub = capacity_bounds)

			if opt_flow.status == 0:
				min_cost = opt_flow.fun
			elif opt_flow.status == 2:
				min_cost = float("inf")
			else:
				raise ValueError("Linear program failed! to find an optimal flow.")

				
		print("Minimum cost ", min_cost, " epsilon feasible ", epsilon_feasible, " epsilon infeasible ", epsilon_infeasible )
		if min_cost == float("inf"):
			raise ValueError("Flow is infeasible.")

		if 1 - min_cost > epsilon_feasible:
			if not use_tree:			
				#optimal_flow_map = recover_optimal_flow(prepared_tree_root)
				#print(optimal_flow_map)
			#else:
				optimal_flow_map = dict([])
				for i in range(len(edges)):
					edge = edges[i]
					optimal_flow_map[edge] = opt_flow.x[i]

					if np.abs(opt_flow.x[i]) not in [0.0,1.0]:
						#print(opt_flow.x)
						raise ValueError("Optimal flow was not integral!!!")

				print(optimal_flow_map)
			if optimal_flow_map != None:
				for e in optimal_flow_map:
					edge = list(e)
					edge.sort()
					#print(edge, optimal_flow_map[e])
			
				coalition_edges = []
				for e in edges:
					if optimal_flow_map[e] != 0:
						coalition_edges.append(e)

				coalition_indicator_vector = np.zeros(len(edges))
				for e in coalition_edges:
					edge_index = edges_to_indices[e]
					coalition_indicator_vector[edge_index] = 1

				coalitions.append(coalition_indicator_vector)

			### Run the linear program to check if the polytope is nonempty.
			linprog_cost = np.ones(len(edges))
			special_edge_indicator = np.zeros(len(edges))
			special_edge_indicator[special_edge_index] = 1
			print("Starting the core linear program. Number of coalitions so far ", len(coalitions))
			if len(coalitions) > max_num_constraints:
				raise ValueError("Exceeded the number of allowed coalition constraints")
			opt_result = optimize.linprog(linprog_cost, A_eq= np.array([np.ones(len(edges)),special_edge_indicator]),
			 b_eq = np.array([1, 0]), A_ub =-np.array(coalitions), b_ub = np.ones(len(coalitions))*(epsilon_feasible-1)   )

			### If the polytope is empty. This is unfeasible. And means we should backtrack 
			if opt_result.status == 2:
				#We should increase epsilon.
				new_epsilon_feasible = 2*(epsilon_feasible-epsilon_infeasible) + epsilon_infeasible
				epsilon_infeasible = epsilon_feasible
				epsilon_feasible = new_epsilon_feasible
				
			## The LP is feasible
			elif opt_result.status == 0:
				imputation = opt_result.x
			else:
				raise ValueError("The linear program failed")

		else:
			### the imputation is a feasible imputation for the epsilon core. We have solved it.
			epsilon_feasible = (epsilon_feasible-epsilon_infeasible)/2.0 + epsilon_infeasible

	return epsilon_feasible

def run_geometric( n, radius, k, use_tree, max_num_constraints = float("inf")):
	G = get_geometric_graph(n=n, radius=radius)
	prunning_cost = 5
	special_edge = [1,len(G.nodes)-1]
	G.add_edge(special_edge[0], special_edge[1])
	special_edge = frozenset(special_edge)
	epsilon = get_epsilon_core( k = k, use_tree = use_tree, 
	special_edge = special_edge, G = G, prunning_cost = prunning_cost, max_num_constraints = max_num_constraints)
	return epsilon


def run_pipelines(small_graph, k, use_tree):
	special_edge = frozenset({'"Poland"', '"Italy"'})
	G = get_graph_pipelines(shortform = small_graph, special_edge = special_edge)
	prunning_cost = 5
	epsilon = get_epsilon_core(k = k, use_tree = use_tree, 
		special_edge = special_edge, G = G, prunning_cost = prunning_cost)
	return epsilon

def run_duplication(n, p, k, use_tree):
	special_edge = (1,19)
	G = get_duplication_graph(n=40, p = .1, special_edge = special_edge )
	special_edge = frozenset(special_edge)
	prunning_cost = 5
	epsilon = get_epsilon_core(k = k, use_tree = use_tree, 
		special_edge = special_edge, G = G, prunning_cost = prunning_cost)
	return epsilon		


def main(r_list = None, n_list = None):
	graph_type  = "geometric"
	k = 1
	use_tree = False

	if graph_type == "pipelines":
		epsilon = run_pipelines(small_graph = False, k = k, use_tree = use_tree)
		print("Epsilon ", epsilon)
	elif graph_type == "geometric":
		trials_per_round =5
		solution_list = []
		for r in r_list:
			for n in n_list:
				for _ in range(trials_per_round):
					print("r ", r, "n ", n)
					try:
						epsilon = run_geometric(n = n, radius = r, k =k, use_tree = use_tree, max_num_constraints = 110)
						solution = ("r", r, "n", n, "epsilon", epsilon )
					except:						
						solution = ("r", r, "n", n, "epsilon", None )
					solution_list.append(solution)
	elif graph_type == "duplication":
		epsilon = run_duplication(n=40, p = .1, k = k, use_tree = use_tree)
	else:
		raise ValueError("Graph type not implemented. ", graph_type)

	return solution_list

if __name__ == "__main__":
	r_list  = [ .6,]
	n_list = [20]

	solution_list = main(r_list = r_list, n_list = n_list)
	def average(n, r):
		count = 0
		av = 0
		for a in solution_list:
			if a[1] == r and a[3] == n and a[5] != None:
				av += a[5]
				if a[5] > -1:
					count += 1
		if count > 0:
			av = av*1.0/count
		else:
			av = None
		return av, count

	averages_list = []
	for r in r_list:
		for n in n_list:
			av, count = average(n,r)
			averages_list.append(("n", n, "r", r, "average", av, "count", count))