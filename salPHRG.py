import os
import sys
import re
import argparse
import traceback
import pandas as pd
import networkx as nx
import tree_decomposition as td
import net_metrics as metrics
import graph_sampler as gs
import probabilistic_cfg as pcfg

__version__ = "0.1.0"
__author__ = ['Salvador Aguinaga', 'Rodrigo Palacios', 'David Chaing', 'Tim Weninger']

# salPHRG does nth inf. mirror

prod_rules = {}
debug = DBG = False


def graph_checks(G):
    ## Target number of nodes
    global num_nodes
    num_nodes = G.number_of_nodes()

    if not nx.is_connected(G):
        print "Graph must be connected";
        os._exit(1)

    if G.number_of_selfloops() > 0:
        print "Graph must be not contain self-loops";
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
    N.append(["S"])  # starting node
    ttt = 0
    # pick non terminal
    num = 0
    for r in rule_list:
        rule = grammar.by_id[r][0]
        lhs_match = matcher(rule.lhs, N)
        e = []
        match = []
        for tup in lhs_match:
            match.append(tup[0])
            e.append(tup[1])
        lhs_str = "(" + ",".join(str(x) for x in sorted(e)) + ")"

        new_idx = {}
        n_rhs = rule.rhs
        # print lhs_str, "->", n_rhs
        for x in n_rhs:
            new_he = []
            he = x.split(":")[0]
            term_symb = x.split(":")[1]
            for y in he.split(","):
                if y.isdigit():  # y is internal node
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
                else:  # y is external node
                    for tup in lhs_match:  # which external node?
                        if tup[1] == y:
                            new_he.append(tup[0])
                            break
            # prod = "(" + ",".join(str(x) for x in new_he) + ")"
            if term_symb == "N":
                N.append(sorted(new_he))
            elif term_symb == "T":
                H.append(new_he)
                #print n_rhs, new_he
        match = sorted(match)
        N.remove(match)

    newG = nx.Graph()
    for e in H:
        if (len(e) == 1):
            newG.add_node(e[0])
        else:
            newG.add_edge(e[0], e[1])

    return newG, D


def probabilistic_hrg(G, num_samples=1):

    graphletG = []

    print G.number_of_nodes()
    print G.number_of_edges()

    G.remove_edges_from(G.selfloop_edges())
    giant_nodes = max(nx.connected_component_subgraphs(G), key=len)
    G = nx.subgraph(G, giant_nodes)

    num_nodes = G.number_of_nodes()

    print G.number_of_nodes()
    print G.number_of_edges()

    graph_checks(G)

    print
    print "--------------------"
    print "-Tree Decomposition-"
    print "--------------------"

    if num_nodes >= 500:
        for Gprime in gs.rwr_sample(G, 2, 100):
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

    print
    print "--------------------"
    print "- Production Rules -"
    print "--------------------"

    for k in prod_rules.iterkeys():
        print k
        s = 0
        for d in prod_rules[k]:
            s += prod_rules[k][d]
        for d in prod_rules[k]:
            prod_rules[k][d] = float(prod_rules[k][d]) / float(s)  # normailization step to create probs not counts.
            #print '\t -> ', d, prod_rules[k][d]

    rules = []
    id = 0
    for k, v in prod_rules.iteritems():
        sid = 0
        for x in prod_rules[k]:
            rhs = re.findall("[^()]+", x)
            rules.append(("r%d.%d" % (id, sid), "%s" % re.findall("[^()]+", k)[0], rhs, prod_rules[k][x]))
            #print ("r%d.%d" % (id, sid), "%s" % re.findall("[^()]+", k)[0], rhs, prod_rules[k][x])
            sid += 1
        id += 1

    g = pcfg.Grammar('S')
    for (id, lhs, rhs, prob) in rules:
        g.add_rule(pcfg.Rule(id, lhs, rhs, prob))

    print "Starting max size"

    g.set_max_size(num_nodes)

    print "Done with max size"

    Hstars = []

    for i in range(0, num_samples):
        rule_list = g.sample(num_nodes)
        # print rule_list
        hstar = grow(rule_list, g)[0]
        # print "H* nodes: " + str(hstar.number_of_nodes())
        # print "H* edges: " + str(hstar.number_of_edges())
        # Hstars.append(hstar)

    return hstar#(Hstars)

def get_phrg_graph_recurrence(graph,recurrence_nbr=1):
  """

  Returns
  -------
  output file of the x,y values
  """
  for r in range(recurrence_nbr):
    g = probabilistic_hrg(graph, 1)
    graph = g

  return graph

def get_hrg_graph_recurrence(graph,recurrence_nbr=1):
  from HRG import stochastic_hrg
  """
  Returns
  -------
  output file of the x,y values
  """
  for r in range(recurrence_nbr):
    g = stochastic_hrg(graph, 1)
    graph = g[0][0]
    print '\t', r, graph.number_of_nodes(),graph.number_of_edges()
  return graph

def get_clgm_graph_recurrence(graph,recurrence_nbr=1):
  """
  Returns
  -------
  output file of the x,y values
  """
  for r in range(recurrence_nbr):
    z = graph.degree().values()
    cl_grph = nx.expected_degree_graph(z)
    graph = cl_grph

  print '\tChung-Lu  -> run:', r, graph.number_of_nodes(),graph.number_of_edges()
  return graph

def grow_graphs_using_krongen(graph, gn, recurrence_nbr=1, graph_vis_bool=False):
  """
  grow graph using krongen given orig graph, gname, and # of recurrences
  Returns
  -------
  nth graph --<kpgm>--
  """
  import math
  from cikm_experiments import kronfit
  from os import environ
  import subprocess

  tsvGraphName = "/tmp/{}kpgraph.tsv".format(gn)
  tmpGraphName = "/tmp/{}kpgraph.tmp".format(gn)
  
  for r in range(recurrence_nbr):
    if environ['HOME'] == '/home/saguinag':
      args = ("/home/saguinag/Software/Snap-2.4/examples/krongen/krongen", "-i:{}".format(tmpGraphName),"-n0:2", "-m:\"0.9 0.6; 0.6 0.1\"", "-gi:5")
    elif environ['HOME'] == '/Users/saguinag':
      args = ("/Users/saguinag/ToolSet/Snap-2.4/examples/krongen/krongen", "-i:{}".format(tmpGraphName),"-n0:2", "-m:\"0.9 0.6; 0.6 0.1\"", "-gi:5")
    else:
      args = ('./kronfit.exe -i:tmp.txt -n0:2 -m:"0.9 0.6; 0.6 0.1" -gi:5')

    kp_graphs = []
    k = int(math.log(graph.number_of_nodes(),2))+1 # Nbr of Iterations
    print 'k:',k,'n',graph.number_of_nodes()
    
    P = kronfit(graph) #[[0.9999,0.661],[0.661,     0.01491]]

    if DBG: print ":recurrence #>",recurrence_nbr

    # KPG = product.kronecker_random_graph(k, P)
    M = '-m:"{} {}; {} {}"'.format(P[0][0], P[0][1], P[1][0], P[1][1])
    if environ['HOME'] == '/home/saguinag':
      args = ("/home/saguinag/Software/Snap-2.4/examples/krongen/krongen", "-o:"+tsvGraphName, M, "-i:{}".format(k))
    elif environ['HOME'] == '/Users/saguinag':
      args = ("/Users/saguinag/ToolSet/Snap-2.4/examples/krongen/krongen", "-o:"+tsvGraphName, M, "-i:{}".format(k))
    else:
      args = ('./krongen.exe -o:{} '.format(tmpGraphName) +M +'-i:{}'.format(k+1))

    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    #output = popen.stdout.read()

    if os.path.exists(tsvGraphName):
      KPG = nx.read_edgelist(tsvGraphName, nodetype=int)
    else:
      print "!! Error, file is missing"

    for u,v in KPG.selfloop_edges():
      KPG.remove_edge(u, v)
    graph = KPG
    if DBG: print 'Avg Deg:', nx.average_degree_connectivity(graph)
    import phoenix.visadjmatrix as vis
    # vis.draw_sns_adjacency_matrix(graph)
    vis.draw_sns_graph(graph)

  # returns the ith recurrence 
  return KPG#kp_graphs

def grow_graphs_using_kpgm(graph,recurrence_nbr=1):
  """
  Returns
  -------
  nth graph --<kpgm>--
  """
  import product
  import math
  from cikm_experiments import kronfit
  kp_graphs = []
  k = int(math.log(graph.number_of_nodes(),2))
  P = kronfit(graph) #[[0.9999,0.661],[0.661,     0.01491]]
  print P
  exit()
  for r in range(recurrence_nbr):
    KPG = product.kronecker_random_graph(k, P)
    for u,v in KPG.selfloop_edges():
      KPG.remove_edge(u, v)
    graph = KPG
    print '\tKronecker -> run:', r, graph.number_of_nodes(),graph.number_of_edges()
    kp_graphs.append(graph)
  return kp_graphs


def get_kpgm_graph_recurrence(graph,recurrence_nbr=1):
  """
  Returns
  -------
  nth graph --<kpgm>--
  """
  import product
  import math
  from cikm_experiments import kronfit
  for r in range(recurrence_nbr):
    k = int(math.log(graph.number_of_nodes(),2))
    P = kronfit(graph) #[[0.9999,0.661],[0.661,     0.01491]]
    KPG = product.kronecker_random_graph(k, P, directed=False)

    for u,v in KPG.selfloop_edges():
      KPG.remove_edge(u, v)
    graph = KPG
    if KPG.is_directed(): print '!!!!!! is directed'
  print '\tKronecker -> run:', r, graph.number_of_nodes(),graph.number_of_edges()
  return KPG


def get_parser():
  parser = argparse.ArgumentParser(description='hrgm: Hyperedge Replacement Grammars Model')
  parser.add_argument('graph', metavar='GRAPH', help='graph path to process')
  parser.add_argument('--version', action='version', version=__version__)
  return parser

def infinity_degree(g, name, nbr_of_averages=1):
  if 0:
    # ---- PHRG ----
    #
    mf = pd.DataFrame()
    for runs in range(nbr_of_averages):
      print '\t>', g.number_of_nodes(), g.number_of_edges()
      synth_g = get_phrg_graph_recurrence(g, 10)

      # After 10 recurrences:
      d = synth_g.degree()
      df = pd.DataFrame.from_dict(d.items())
      gb = df.groupby([1]).count()
      mf = pd.concat([mf, gb], axis=1)

    mf['pk'] = mf.mean(axis=1)/float(g.number_of_nodes())
    mf['k'] = mf.index.values
    #print mf
    out_tsv = '../Results/'+name+"_phrg_degree.tsv"
    with open(out_tsv, "w") as f: f.write('k\tpk\n')
    mf[['k','pk']].to_csv(out_tsv, sep='\t', index=False, header=False, mode="w")

  #
  # ---- CHUNG LU ----
  #
  mf = pd.DataFrame()
  for run in range(nbr_of_averages):
    print '\t<',run,'>', g.number_of_nodes(), g.number_of_edges()
    synth_g = get_clgm_graph_recurrence(g, 10)

    # After 10 recurrences:
    d = synth_g.degree()
    df = pd.DataFrame.from_dict(d.items())
    gb = df.groupby([1]).count()
    mf = pd.concat([mf, gb], axis=1)

  mf['pk'] = mf.mean(axis=1)/float(g.number_of_nodes())
  mf['k'] = mf.index.values
  #print mf
  out_tsv = '../Results/'+name+"_clgm_degree.tsv"
  mf[['k','pk']].to_csv(out_tsv, sep='\t', index=False, header=True, mode="w")

  #
  # ---- Kronecker ----
  #
  mf = pd.DataFrame()
  for run in range(nbr_of_averages):
    print '\t<',run,'>', g.number_of_nodes(), g.number_of_edges()
    synth_g = grow_graphs_using_krongen(g, 10)

    # After 10 recurrences:
    d = synth_g.degree()
    df = pd.DataFrame.from_dict(d.items())
    gb = df.groupby([1]).count()
    mf = pd.concat([mf, gb], axis=1)

  mf['pk'] = mf.mean(axis=1)/float(g.number_of_nodes())
  mf['k'] = mf.index.values

  out_tsv = '../Results/'+name+"_kpgm_degree.tsv"
  mf[['k','pk']].to_csv(out_tsv, sep='\t', index=False, header=True, mode="w")

def infinity_clustering_coefficients(g, name, nbr_of_averages=1):
  # ---- ORIG ----
  #
  df = pd.DataFrame.from_dict(g.degree().items())
  df.columns=['v','k']
  cf = pd.DataFrame.from_dict(nx.clustering(g).items())
  cf.columns=['v','cc']
  df = pd.merge(df,cf,on='v')
  gb = df.groupby(['k']).count()
  gb['avg_cc'] =df.groupby(['k']).mean()
  out_tsv = '../Results/'+name+"_orig_clustering.tsv"
  gb[['avg_cc']].to_csv(out_tsv, sep="\t", header=True, index=True)

  if 0:
    # ---- PHRG ----
    mf = pd.DataFrame()
    for run in range(nbr_of_averages):
      print '\t<',run,'>', g.number_of_nodes(), g.number_of_edges()
      synth_g = get_phrg_graph_recurrence(g, 2)

      # After 10 recurrences
      df = pd.DataFrame.from_dict(synth_g.degree().items())
      df.columns=['v','k']
      cf = pd.DataFrame.from_dict(nx.clustering(synth_g).items())
      cf.columns=['v','cc']
      df = pd.merge(df,cf,on='v')
      # gb = df.groupby(['k']).count()
      # gb['avg_cc'] =df.groupby(['k']).mean()
      mf = pd.concat([mf, df])
    gb = mf.groupby(['k']).mean()

    out_tsv = '../Results/'+name+"_phrg_clustering.tsv"
    gb[['cc']].to_csv(out_tsv, sep="\t", header=True, index=True)



  # ---- CHUNG LU ----
  #
  mf = pd.DataFrame()
  for run in range(nbr_of_averages):
    print '\t<',run,'>', g.number_of_nodes(), g.number_of_edges()
    synth_g = get_clgm_graph_recurrence(g, 10)

    # After 10 recurrences:
    df = pd.DataFrame.from_dict(synth_g.degree().items())
    df.columns=['v','k']
    cf = pd.DataFrame.from_dict(nx.clustering(synth_g).items())
    cf.columns=['v','cc']
    df = pd.merge(df,cf,on='v')
    # gb = df.groupby(['k']).count()
    # gb['avg_cc'] =df.groupby(['k']).mean()
    mf = pd.concat([mf, df])
  gb = mf.groupby(['k']).mean()

  out_tsv = '../Results/'+name+"_clgm_clustering.tsv"
  gb[['cc']].to_csv(out_tsv, sep="\t", header=True, index=True)

  #
  # ---- Kronecker ----
  #
  mf = pd.DataFrame()
  for run in range(nbr_of_averages):
    print '\t<',run,'>', g.number_of_nodes(), g.number_of_edges()
    synth_g = get_kpgm_graph_recurrence(g, 10)

    # After 10 recurrences:
    df = pd.DataFrame.from_dict(synth_g.degree().items())
    df.columns=['v','k']
    #print type(synth_g), synth_g.number_of_nodes(), synth_g.number_of_edges()
    cf = pd.DataFrame.from_dict(nx.clustering(synth_g).items())
    cf.columns=['v','cc']
    df = pd.merge(df,cf,on='v')
    # gb = df.groupby(['k']).count()
    # gb['avg_cc'] =df.groupby(['k']).mean()
    mf = pd.concat([mf, df])
    # print df.groupby(['k']).mean()
  gb = mf.groupby(['k']).mean()

  out_tsv = '../Results/'+name+"_kpgm_clustering.tsv"
  gb[['cc']].to_csv(out_tsv, sep="\t", header=True, index=True)


def basic_network_statistics(G, name=""):
  import pprint 
  # print '\degree_probability_distribution','-'*40
  # metrics.draw_degree_probability_distribution(orig_g_M=[cg], HRG_M=[], pHRG_M=[], chunglu_M=[], kron_M=[])
  # print '\tdraw_network_value','-'*40
  # metrics.draw_network_value(orig_g_M=[cg], HRG_M=[], pHRG_M=[], chunglu_M=[], kron_M=[])
  # print '\tdraw_hop_plot','-'*40
  # metrics.draw_hop_plot(orig_g_M=[cg], HRG_M=[], pHRG_M=[], chunglu_M=[], kron_M=[])
  # print '\tdraw_kcore_decomposition','-'*40
  # metrics.draw_kcore_decomposition(orig_g_M=[cg], HRG_M=[], pHRG_M=[], chunglu_M=[], kron_M=[])
  # print '\tdraw_clustering_coefficients','-'*40
  # metrics.draw_clustering_coefficients(orig_g_M=[cg], HRG_M=[], pHRG_M=[], chunglu_M=[], kron_M=[])
  # print '\tdraw_assortativity_coefficients','-'*40
  # metrics.draw_assortativity_coefficients(orig_g_M=[cg], HRG_M=[], pHRG_M=[], chunglu_M=[], kron_M=[])
  ofname = '../Results/basic_stats_{}.log'.format(name)
  stats = {}
  with open(ofname, 'w') as f:
    
    
    # deg dist
    d = G.degree()
    df = pd.DataFrame.from_dict(d.items())
    stats['meandeg'] = df[1].mean()
    stats['maxdeg'] = df[1].max()
    stats['avgcc'] = nx.average_clustering(G)
 
  pprint.pprint(stats) 
  return

def basic_statistics_on_generated_networks(G, name="", nbr_of_averages=1):
  # ---- PHRG ----
  #
  mf = pd.DataFrame()
  for runs in range(nbr_of_averages):
    synth_g = get_phrg_graph_recurrence(G, 1)
    df = pd.DataFrame.from_dict(synth_g.degree().items())
    df.columns = ['v','k']
    mf = pd.concat([mf,df])
  print "basic_statistics_on_generated_networks:", name, mf['k'].mean()

def compute_infinity_mirror_graph_sets(cg, name):
  import shelve

  ithrun = 10
  for r in range(1, ithrun + 1):
    cl_graphs = []
    kp_graphs = []

    for n in range(10):
      print '~~ run:',n,'(avg)','~'*20
      cl = get_clgm_graph_recurrence(cg, r)
      cl_graphs.append(cl)
      kp = grow_graphs_using_krongen(cg, name, r) # get the rth synth graph
      kp_graphs.append(kp)

    print len(cl_graphs),'::',len(kp_graphs)

    shl = shelve.open("../Results/synthg_{}_10x{}th.shl".format(name, r)) # the same filename tha    t you used before, please
    shl['clgm'] = (cl_graphs, name)
    shl['kpgm'] = (kp_graphs, name)
    shl.close()




def main():
  global name

  parser = get_parser()
  args = vars(parser.parse_args())

  if not args['graph']:
    parser.print_help()
    os._exit(1)
  
  try:
    cg = nx.read_edgelist(args['graph'])
  except Exception, e:
    print 'ERROR, UNEXPECTED READ EXCEPTION'
    print str(e)
    try: 
      g = nx.read_edgelist(args['graph'], comments="%")
    except Exception, e:
      print 'Error, unexpected read edgelist exception'
      print str(e)
      
  print 'Read Edgelist File'
  name = os.path.basename(args['graph']).rstrip('.txt')

  # infinity_degree(cg,name, nbr_of_averages=10)
  # infinity_clustering_coefficients(cg, name, nbr_of_averages=10)
  #basic_network_statistics(cg, name)
  #basic_statistics_on_generated_networks(cg, name, 4)
  compute_infinity_mirror_graph_sets(cg, name)
  # basic_network_statistics(cg, name)
  # basic_statistics_on_generated_networks(cg, name, 4)
  # basic_network_statistics(cg, name)
  # basic_statistics_on_generated_networks(cg, name, 4)
  #basic_statistics_on_generated_networks(cg, name, 4)


if __name__ == "__main__":
  try:
    main()
  except Exception, e:
    print 'ERROR, UNEXPECTED SAVE PLOT EXCEPTION'
    print str(e)
    traceback.print_exc()
    os._exit(1)
  sys.exit(0)
