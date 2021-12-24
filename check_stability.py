import re
import math
import ppl
import gmpy2
import ppcd_input
import copy
import networkx as nx
import matplotlib.pyplot as plt
import ppcd_abstract
import numpy as np

def declare_vars(var_list):
    for i in range(len(var_list)):
        define_str = "global " + var_list[i] + "\n" +\
             var_list[i] + "=" + "ppl.Variable(" + str(i) + ")"
        #print(define_str)
        exec(define_str)

#check if there is a negative cycle in graph from source using Bellman-Ford
#if yes return the cycle otherwise return False
def isNegCycle(graph,source):
    dist = {}
    path = {}
    V = list(graph.nodes())
    E = list(graph.edges.data())
    for v in V:
        dist[v] = math.inf
        path[v] = []
    dist[source] = 0
    path[source] = [source]
    for i in range(1,len(V)):
        for j in range(len(E)):
            u = E[j][0]
            v = E[j][1]
            weight = E[j][2]['weight']
            if dist[u] != math.inf and dist[u] + weight < dist[v]:
                dist[v] = dist[u] + weight
                path[v] = path[u] + [v]

    #print(dist)
    #print(path)
    for i in range(len(E)):
        u = E[i][0]
        v = E[i][1]
        weight = E[i][2]['weight']
        if dist[u] != math.inf and dist[u] + weight < dist[v]:
            path[v] = path[u] + [v]
            for j in range(len(path[v])):
                if path[v][j] in path[v][:j]: break
            #print(j)
            return(path[v][:(j+1)])
        
    return False

'''
G = nx.DiGraph()
G.add_edge(0,1,weight=2)
G.add_edge(1,2,weight=-2)
G.add_edge(2,0,weight=0)
G.add_edge(3,0,weight=1.5)
print(isNegCycle(G,0))
'''

#Check for infinite weight edge or cycle with with weight>1 in G from source
def check_abs_stability(G,source):
    G_anal = ppcd_abstract.graph_for_analysis(G,source)
    #print(G_anal.edges())
    for e in list(G_anal.edges.data()):
        if e[2]['weight'] == None:
            #print("Exploding region reachable from source.")#
            return -1

    truth_value = isNegCycle(G_anal,source)
    if truth_value == False:
        #print("Abstract system stable.")#
        return 1
    else:
        #print("Abstract system not stable.")#
        return 0

'''
PPCD = ppcd_input.readPPCD("ppcd_example1.txt")
G,facet_list = ppcd_abstract.create_abs_graph(PPCD)
source = list(G.nodes())[0]
#G_anal = ppcd_abstract.graph_for_analysis(G, source)
#ppcd_abstract.draw_Ganal(G_anal, facet_list)
check_abs_stability(G,source)
'''

def find_stationary_dist(G):
    var_list = []
    V = list(G.nodes())
    E = list(G.edges.data())
    for i in range(len(V)):
        var_list.append("x"+str(i))
    
    declare_vars(var_list)

    stationary_eqs = ppl.Constraint_System()
    sum_expr = {"const" : 1}
    for j in range(len(V)):
        expr = {"const" : 0}
        denom_lcm = 1
        for e in E:
            if e[1] == V[j]:
                if e[0] == V[j]: 
                    expr[var_list[j]] = (e[2]['prob'] - 1)
                    denom = gmpy2.denom(expr[var_list[j]])
                    denom_lcm = math.lcm(denom_lcm,denom)
                else:
                    i = V.index(e[0])
                    expr[var_list[i]] = e[2]['prob']
                    denom = gmpy2.denom(expr[var_list[i]])
                    denom_lcm = math.lcm(denom_lcm,denom)
        #print(denom_lcm)
        for key in expr: 
            new_coeff = expr[key]*denom_lcm
            #print(new_coeff)
            expr[key] = new_coeff

        ineq = ppcd_abstract.build_ineq_from_expr(expr, "==")
        #print(eval(ineq))
        sum_expr[var_list[j]] = 1
        stationary_eqs.insert(eval(ineq))

    ineq = ppcd_abstract.build_ineq_from_expr(sum_expr, "==")
    #print(ineq)
    stationary_eqs.insert(eval(ineq))
    #print(stationary_eqs)
    
    stationary_poly = ppl.C_Polyhedron(stationary_eqs)
    #print(stationary_poly)
    
    stationary_dist = {}
    if stationary_poly.is_empty() == False:
        for i in range(len(V)):
            expr_to_max = "1*"+var_list[i]
            target = stationary_poly.maximize(eval(expr_to_max))
            if target['bounded']: stationary_dist[V[i]] = \
                gmpy2.mpq(target['sup_n'],target['sup_d'])
    #print(stationary_dist)
    return stationary_dist


'''
PPCD = ppcd_input.readPPCD("ppcd_example1.txt")
G,facet_list = ppcd_abstract.create_abs_graph(PPCD)
source = list(G.nodes())[0]
G_anal = ppcd_abstract.graph_for_analysis(G, source)
#print(G_anal.edges.data())
stationary_dist = find_stationary_dist(G_anal)
'''


def check_as_stability(G,source):
    G_anal = ppcd_abstract.graph_for_analysis(G,source)
    for e in list(G_anal.edges.data()):
        if e[2]['weight'] == None:
            #print("System is not almost surely stable")
            return 0

    stationary_dist = find_stationary_dist(G_anal)
    as_avg_weight = 0
    for e in list(G_anal.edges.data()):
        as_avg_weight = as_avg_weight - stationary_dist[e[0]]*(e[2]['weight'])

    #print(as_avg_weight)
    
    if as_avg_weight <= 0 :
        #print("System is almost surely stable")
        return 1
    else:
        #print("System is not almost surely stable")
        return 0

'''
PPCD = ppcd_input.readPPCD("ppcd_example3.txt")
G,facet_list = ppcd_abstract.create_abs_graph(PPCD)
source = list(G.nodes())[0]
check_as_stability(G, source)
'''