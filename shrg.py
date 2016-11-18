# class $class $($object$):
#   """$cls_doc$"""
#
#   def __init__(self, $args$):
#     """Constructor for $class$"""
#     $END$
import networkx as nx
import pandas as pd
import re

_author_ = 'saguinag'

class HRG:
  def __init__(self, edglst_fname):
    """Constructor for $class$"""
    self.edglst_fname = edglst_fname

  def learngrammar(self):
    g = nx.read_edgelist(self.edglst_fname,
                         comments="%")
    print nx.info(g)

class GrammarEM:
  def __init__(self, fname):
    """Constructor for $class$"""
    self.fname = fname

  def reformat_prod_rules(self,in_grammar_lst):
    # for nbr,pr,p in in_grammar_lst.__iter__():
    #   print nbr,pr,p
    #   break
    for p in in_grammar_lst:
      print len(p)
      break

  def load_grammar(self):
    df = pd.read_csv(self.fname,sep="\t", header=None)
    df[1]= df[1].apply(lambda x: x.translate(None, "(){}"))
    df['lhs'] = df[1].apply(lambda x: x.split()[0].rstrip(',') )
    df['rhs'] = df[1].apply(lambda x: list(x.split('[', 1)[1].split(']')[0].split(',')) )
    df = df[[0,'lhs','rhs',2]]
    print df.head()

#HRG('/Users/saguinag/Research/Phoenix/demo_graphs/subelj_euroroad/out.subelj_euroroad_euroroad').learngrammar()

GrammarEM('outfile.tsv').load_grammar()