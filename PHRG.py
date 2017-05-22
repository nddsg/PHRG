import os
import re
import networkx as nx
import probabilistic_cfg as pcfg
import graph_sampler as gs
import tree_decomposition as td
import random
import probabilistic_gen as pg

# prod_rules = {}
DEBUG = False


def graph_checks(G):
		## Target number of nodes
		global num_nodes
		num_nodes = G.number_of_nodes()

		if not nx.is_connected(G):
				if DEBUG: print "Graph must be connected";
				os._exit(1)

		if G.number_of_selfloops() > 0:
				if DEBUG: print "Graph must be not contain self-loops";
				os._exit(1)


def matcher(lhs, N):
		if lhs == "S":
				return [("S", "S")]
		m = []
		for x in N:
				if len(x) == lhs.count(",") + 1:
						i = 0
						for y in lhs.split(","):
								m.append((x[i], y))
								i += 1
						return m


def binarize(tree):
		(node, children) = tree
		children = [binarize(child) for child in children]
		if len(children) <= 2:
				return (node, children)
		else:
				# Just make copies of node.
				# This is the simplest way to do it, but it might be better to trim unnecessary members from each copy.
				# The order that the children is visited is arbitrary.
				binarized = (node, children[:2])
				for child in children[2:]:
						binarized = (node, [binarized, child])
				return binarized


def grow(rule_list, grammar, diam=0):
		D = list()
		newD = diam
		H = list()
		N = list()
		N.append(["S"])	# starting node
		ttt = 0
		# pick non terminal
		num = 0
		for r in rule_list:
				rule = grammar.by_id[r][0]
				lhs_match = matcher(rule.lhs, N)
				e = []	# edge list
				match = []
				for tup in lhs_match:
						match.append(tup[0])
						e.append(tup[1])
				lhs_str = "(" + ",".join(str(x) for x in sorted(e)) + ")"

				new_idx = {}
				n_rhs = rule.rhs
				if 0: print lhs_str, "->", n_rhs
				for x in n_rhs:
						new_he = []
						he = x.split(":")[0]
						term_symb = x.split(":")[1]
						for y in he.split(","):
								if y.isdigit():	# y is internal node
										if y not in new_idx:
												new_idx[y] = num
												num += 1
												if diam > 0 and num >= newD and len(H) > 0:
														newD = newD + diam
														newG = nx.Graph()
														for e in H:
																if (len(e) == 1):
																		newG.add_node(e[0])
																else:
																		newG.add_edge(e[0], e[1])
																		# D.append(metrics.bfs_eff_diam(newG, 20, 0.9))
										new_he.append(new_idx[y])
								else:	# y is external node
										for tup in lhs_match:	# which external node?
												if tup[1] == y:
														new_he.append(tup[0])
														break
						# prod = "(" + ",".join(str(x) for x in new_he) + ")"
						if term_symb == "N":
								N.append(sorted(new_he))
						elif term_symb == "T":
								H.append(new_he)	# new (h)yper(e)dge
								# print n_rhs, new_he
				match = sorted(match)
				N.remove(match)

		newG = nx.Graph()
		for e in H:
				if (len(e) == 1):
						newG.add_node(e[0])
				else:
						newG.add_edge(e[0], e[1])

		return newG, D


def probabilistic_hrg(G, n=None):
		'''
			Rule extraction procedure

					'''
		if G is None: return

		G.remove_edges_from(G.selfloop_edges())
		giant_nodes = max(nx.connected_component_subgraphs(G), key=len)
		G = nx.subgraph(G, giant_nodes)

		if n is None:
				num_nodes = G.number_of_nodes()
		else:
				num_nodes = n

		graph_checks(G)

		if DEBUG: print
		if DEBUG: print "--------------------"
		if DEBUG: print "-Tree Decomposition-"
		if DEBUG: print "--------------------"
		prod_rules = {}
		if num_nodes >= 500:
				for Gprime in gs.rwr_sample(G, 2, 300):
						T = td.quickbb(Gprime)
						root = list(T)[0]
						T = td.make_rooted(T, root)
						T = binarize(T)
						root = list(T)[0]
						root, children = T
						td.new_visit(T, G, prod_rules)
		else:
				T = td.quickbb(G)
				root = list(T)[0]
				T = td.make_rooted(T, root)
				T = binarize(T)
				root = list(T)[0]
				root, children = T
				td.new_visit(T, G, prod_rules)

		if DEBUG: print
		if DEBUG: print "--------------------"
		if DEBUG: print "- Production Rules -"
		if DEBUG: print "--------------------"

		for k in prod_rules.iterkeys():
				if DEBUG: print k
				s = 0
				for d in prod_rules[k]:
						s += prod_rules[k][d]
				for d in prod_rules[k]:
						prod_rules[k][d] = float(prod_rules[k][d]) / float(s)	# normailization step to create probs not counts.
						if DEBUG: print '\t -> ', d, prod_rules[k][d]

		# pp.pprint(prod_rules)

		rules = []
		id = 0
		for k, v in prod_rules.iteritems():
				sid = 0
				for x in prod_rules[k]:
						rhs = re.findall("[^()]+", x)
						rules.append(("r%d.%d" % (id, sid), "%s" % re.findall("[^()]+", k)[0], rhs, prod_rules[k][x]))
						if DEBUG: print ("r%d.%d" % (id, sid), "%s" % re.findall("[^()]+", k)[0], rhs, prod_rules[k][x])
						sid += 1
				id += 1

		return rules

def probabilistic_hrg_learning(G, num_samples=1, n=None, prod_rules=None):

		graphletG = []

		# print G.number_of_nodes()
		# print G.number_of_edges()

		G.remove_edges_from(G.selfloop_edges())
		giant_nodes = max(nx.connected_component_subgraphs(G), key=len)
		G = nx.subgraph(G, giant_nodes)

		if n is None:
				num_nodes = G.number_of_nodes
		else:
				num_nodes = n

		# print G.number_of_nodes()
		# print G.number_of_edges()

		graph_checks(G)

		# print
		# print "--------------------"
		# print "-Tree Decomposition-"
		# print "--------------------"

		if num_nodes >= 500:
				for Gprime in gs.rwr_sample(G, 2, 300):
						T = td.quickbb(Gprime)
						root = list(T)[0]
						T = td.make_rooted(T, root)
						T = binarize(T)
						root = list(T)[0]
						root, children = T
						td.new_visit(T, G, prod_rules)
		else:
				T = td.quickbb(G)
				root = list(T)[0]
				T = td.make_rooted(T, root)
				T = binarize(T)
				root = list(T)[0]
				root, children = T
				td.new_visit(T, G, prod_rules)
		# print 'root', [x for x in T[0]]#, type(root)
		# import pprint as pp
		# pp.pprint([x for x in T])
		'''
		for x in T:
			if isinstance(x,(frozenset)):
				print '\t',x
			else:
				print [type(s) for s in x if isinstance(x,(list))]
		'''
		##while isinstance(T,(tuple,list,)) and len(T):
		##	for x in T:
		##		if isinstance(x,(frozenset)):
		##			print'\t',	x
		##		else:
		##			T = x


		# print
		# print "--------------------"
		# print "- Production Rules -"
		# print "--------------------"

		for k in prod_rules.iterkeys():
				# print k
				s = 0
				for d in prod_rules[k]:
						s += prod_rules[k][d]
				for d in prod_rules[k]:
						prod_rules[k][d] = float(prod_rules[k][d]) / float(s)	# normailization step to create probs not counts.
						# print '\t -> ', d, prod_rules[k][d]

		rules = []
		id = 0
		for k, v in prod_rules.iteritems():
				sid = 0
				for x in prod_rules[k]:
						rhs = re.findall("[^()]+", x)
						rules.append(("r%d.%d" % (id, sid), "%s" % re.findall("[^()]+", k)[0], rhs, prod_rules[k][x]))
						# print ("r%d.%d" % (id, sid), "%s" % re.findall("[^()]+", k)[0], rhs, prod_rules[k][x])
						sid += 1
				id += 1

		return rules

def try_combination(lhs, N):
		for i in range(0, len(N)):
				n = N[i]
				if lhs[0] == "S":
						break
				if len(lhs) == len(n):
						t = []
						for i in n:
								t.append(i)
						random.shuffle(t)
						return zip(t, lhs)
		return []

if __name__ == "__main__":

		# Example From PAMI Paper
		# Graph is undirected
		G = nx.Graph()
		G.add_edge(1, 2)
		G.add_edge(2, 3)
		G.add_edge(2, 4)
		G.add_edge(3, 4)
		G.add_edge(3, 5)
		G.add_edge(4, 6)
		G.add_edge(5, 6)
		G.add_edge(1, 5)

		num_nodes = G.number_of_nodes()


		prod_rules = {}
		p_rules = probabilistic_hrg(G)

		g = pcfg.Grammar('S')
		for (id, lhs, rhs, prob) in p_rules:
				g.add_rule(pcfg.Rule(id, lhs, rhs, prob, True))
		print '> prod rules added to Grammar g'	#
		g.set_max_size(num_nodes*4)
		print '> max-size set.'

		rids = g.sample(num_nodes*4)
		print rids

		new_graph = pg.gen(rids, g)

		print "nodes: " , new_graph.number_of_nodes()
		print "edges: " , new_graph.number_of_edges()
		print
		print



