import random
import networkx as nx


def try_combination(lhs, N):
    lhs = lhs.split(',')
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


n_distribution = {}


def gen(rids, g):
    N = list()
    H = list()
    num = 0
    for r in rids:
        applied_rule = r
        lhs = g.by_id[r][0].lhs
        rhs = g.by_id[r][0].rhs

        lhs_match = try_combination(lhs, N)

        match = []
        for tup in lhs_match:
            match.append(tup[0])

        new_idx = {}
        # print lhs, "->", rhs
        for x in rhs:
            new_he = []
            he = x.split(":")[0]
            term_symb = x.split(":")[1]
            for y in he.split(","):
                if y.isdigit():  # y is internal node
                    if y not in new_idx:
                        new_idx[y] = num
                        num += 1
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

        if len(match) > 0: N.remove(sorted(match))

    newG = nx.Graph()
    for e in H:
        if (len(e) == 1):
            newG.add_node(e[0])
        else:
            newG.add_edge(e[0], e[1])

    return newG
