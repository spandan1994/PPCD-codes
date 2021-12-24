import re
import math
import ppl
import gmpy2
import ppcd_input
import copy
import networkx as nx
import matplotlib.pyplot as plt

'''
Read a linear expression and return as dictionary
Example : x+y+3 will be returned as 
{'const': 3.0, 'x': 1.0, 'y': 1.0}
'''
def readExpr(equation):
    d = {"const":0.0}
    tempeqn = equation
    pattern1 = "\s*([0-9+-.]*)([a-zA-Z]+[0-9_']*)\s*" #identify var
    pattern2 = "\s*([0-9+-.\s]*)\s*" #identify const
    while tempeqn != "":
        try:
            result = re.search(pattern1,tempeqn)
            if(result.group(1) == "" or result.group(1) == "+"):
                value = 1.0
            elif(result.group(1) == "-"):
                value = -1.0
            else:
                value = (eval(result.group(1)))*1.0

            d[result.group(2)] = value
            #print(result.group(2))
            #print(result.group(1))
            tempeqn = re.sub(pattern1,"",tempeqn,1) #replace identified const with null str
        except:
            result = re.search(pattern2,tempeqn)
            #print(result.group(1))
            d["const"] = (eval(result.group(1)))*1.0
            tempeqn = re.sub(pattern2,"",tempeqn,1) #replace identified pattern with null string
    
    return d

#s = " x +y -  3"
#print(readExpr(s))


'''
Read a file with a separate linear expression in each line
and return a list of corresponding dictionaries
'''
def readExprFile(filename):
    fp = open(filename,'r')
    fplines = fp.readlines()

    expr_list = [] #stores linear exprs as dictionaries
    for line in fplines:
        modline = re.sub("\s","",line)
        if modline != '':
            d = readExpr(modline)
            expr_list.append(d)

    return expr_list

#l = readExprFile("equations.txt")
#print(l)

#declare the variables to be used in the constraint system
def declare_vars(var_list):
    for i in range(len(var_list)):
        define_str = "global " + var_list[i] + "\n" +\
             var_list[i] + "=" + "ppl.Variable(" + str(i) + ")"
        #print(define_str)
        exec(define_str)

#builds an inequality from an expr stored as dict
def build_ineq_from_expr(expr,sgn):
    ineq = ""
    for key in expr:
        if key != "const":
            ineq = ineq + str(expr[key]) + "*" + key + "+"
    
    ineq = ineq[:-1] + sgn + str(expr["const"])
    return(ineq)

'''
expr = {'const': 0.0, 'x': 1.0, 'y': 1.0}
print(build_ineq_from_expr(expr,"<="))
>>1.0*x+1.0*y<=0.0
'''

#checks whether a constraint is compatible with a constraint system
def check_compatibility(cs,ineq):
    temp_cs = copy.deepcopy(cs)
    temp_cs.insert(ineq)
    #print(cs); print(temp_cs)

    p = ppl.NNC_Polyhedron(temp_cs)
    #print(p.is_empty())
    if p.is_empty(): return False
    else: return True

'''
x = ppl.Variable(0)
y = ppl.Variable(1)
cs = ppl.Constraint_System()
cs.insert(x<=0)
cs.insert(x+y>=0)
ineq = eval("y<0")
print(check_compatibility(cs,ineq))
#>>False
'''

#Take a non-empty list of Constraint_System()
#Update list with constraints formed by exprs in expr_list[m:n-1]
def update_Constraint_list(cs_list,expr_list,m,n):
    for i in range(m,n):
        new_cs_list = []
        for j in range(len(cs_list)):
            cs1 = cs_list.pop()
            #print(cs1)
            cs2 = copy.deepcopy(cs1)
            ineq1 = eval(build_ineq_from_expr(expr_list[i],">"))
            #print(ineq1)
            ineq2 = eval(build_ineq_from_expr(expr_list[i],"<"))
            #print(ineq2)

            if check_compatibility(cs1,ineq1):
                cs1.insert(eval(build_ineq_from_expr(expr_list[i],">=")))
                new_cs_list.append(cs1)
                #print(new_partition_list)
            if check_compatibility(cs2,ineq2):
                cs2.insert(eval(build_ineq_from_expr(expr_list[i],"<=")))
                new_cs_list.append(cs2)
                #print(new_partition_list)
                
        cs_list = new_cs_list
        #print(partition_list)
    
    return cs_list



#create partitions from a set of lin exprs using ppl and returns list
#of partitions (Constraint_Systems())
def createPartsPPL(expr_list,var_list):
    declare_vars(var_list)

    partition_list = []

    cs_upper = ppl.Constraint_System()
    ineq = build_ineq_from_expr(expr_list[0],">=")
    cs_upper.insert(eval(ineq))
    partition_list.append(cs_upper)

    cs_lower = ppl.Constraint_System()
    ineq = build_ineq_from_expr(expr_list[0],"<=")
    cs_lower.insert(eval(ineq))
    partition_list.append(cs_lower)
    #print(partition_list)
    
    partition_list = update_Constraint_list(partition_list,expr_list,1,len(expr_list))

    return partition_list

'''
expr_list = readExprFile("equations1.txt")
#print(expr_list)
l = createPartsPPL(expr_list,['x','y'])
print(l)
#>>[Constraint_System {x0>=0, -x1>=0, x0+x1>=0}, Constraint_System {x0>=0, -x1>=0, -x0-x1>=0}, 
#Constraint_System {x0>=0, x1>=0, x0+x1>=0}, Constraint_System {-x0>=0, -x1>=0, -x0-x1>=0}, 
#Constraint_System {-x0>=0, x1>=0, x0+x1>=0}, Constraint_System {-x0>=0, x1>=0, -x0-x1>=0}]
'''

#create facets for equality in i-th expr
def createFacetPPL(expr_list,i,var_list):
    declare_vars(var_list)
    
    temp_expr_list = copy.deepcopy(expr_list)
    facet_list = []

    cs = ppl.Constraint_System()
    ineq = build_ineq_from_expr(temp_expr_list.pop(i),"==")
    cs.insert(eval(ineq))
    facet_list.append(cs)

    facet_list = update_Constraint_list(facet_list,temp_expr_list,0,len(temp_expr_list))

    return facet_list

'''
expr_list = readExprFile("equations1.txt")
#print(expr_list)
l = createFacetPPL(expr_list,2,['x','y'])
print(l)
#>>[Constraint_System {x0+x1==0, -x0>=0, x1>=0}, 
#Constraint_System {x0+x1==0, x0>=0, -x1>=0}]
'''

#create all facets
def createAllFacetsPPL(expr_list,var_list):
    declare_vars(var_list)

    facets = []
    for i in range(len(expr_list)):
        facets_i = createFacetPPL(expr_list,i,var_list)
        facets = facets + facets_i

    return facets

'''
expr_list = readExprFile("equations1.txt")
#print(expr_list)
l = createAllFacetsPPL(expr_list,['x','y'])
print(l)
#>>[Constraint_System {x0==0, -x1>=0, -x0-x1>=0}, Constraint_System {x0==0, x1>=0, x0+x1>=0}, 
#Constraint_System {x1==0, -x0>=0, -x0-x1>=0}, Constraint_System {x1==0, x0>=0, x0+x1>=0}, 
#Constraint_System {x0+x1==0, -x0>=0, x1>=0}, Constraint_System {x0+x1==0, x0>=0, -x1>=0}]
'''

#Create constraint system for flow
def create_constraints_flow(flow_constraint_list,var_list,new_var_list):
    flow_cs = ppl.Constraint_System()
    for constraint_dict in flow_constraint_list:
        new_constraint_dict = {"const":0}
        for key in constraint_dict:
            if key in var_list:
                new_key = "(" + new_var_list[var_list.index(key)+len(var_list)] +\
                    "-" + new_var_list[var_list.index(key)] + ")"
                new_constraint_dict[new_key] = constraint_dict[key]
            elif key == "const":
                new_constraint_dict["t"] = -(constraint_dict[key])
            else:
                if constraint_dict[key] == "=":
                    new_sgn = "=="
                else:
                    new_sgn = constraint_dict[key]
        ineq = build_ineq_from_expr(new_constraint_dict,new_sgn)
        #print(ineq)
        flow_cs.insert(eval(ineq))
    #print(flow_cs)

    return flow_cs


def create_constraints_inv(inv_constraint_list,var_list,new_var_list):
    inv_cs = ppl.Constraint_System()
    for constraint_dict in inv_constraint_list:
        new_constraint_dict1 = {"const":0}
        new_constraint_dict2 = {"const":0}
        for key in constraint_dict:
            if key in var_list:
                new_key1 = new_var_list[var_list.index(key)+len(var_list)]
                new_key2 = new_var_list[var_list.index(key)]
                new_constraint_dict1[new_key1] = constraint_dict[key]
                new_constraint_dict2[new_key2] = constraint_dict[key]
            elif key == "const":
                new_constraint_dict1["t"] = -(constraint_dict[key])
                new_constraint_dict2["t"] = -(constraint_dict[key])
            else:
                if constraint_dict[key] == "=":
                    new_sgn = "=="
                else:
                    new_sgn = constraint_dict[key]
        ineq1 = build_ineq_from_expr(new_constraint_dict1,new_sgn)
        ineq2 = build_ineq_from_expr(new_constraint_dict2,new_sgn)
        #print(ineq1)
        inv_cs.insert(eval(ineq1)); inv_cs.insert(eval(ineq2))
    #print(inv_cs)

    return inv_cs


def check_explosion(cs1,additional_cs,expr,flow_constraint_list,inv_constraint_list,var_list):
    new_var_list = []
    for var in var_list:
        new_var_list.append(var+"1")
    for var in var_list:
        new_var_list.append(var+"2")
    new_var_list.append("t")
    #print(new_var_list)
    declare_vars(new_var_list)

    #Create polyhedron for cs1
    p1 = ppl.C_Polyhedron(len(var_list),'universe')#
    p1.add_constraints(cs1)

    #Create polyhedron for cs2
    p2 = ppl.C_Polyhedron(len(var_list),'universe')#
    #p2.add_constraints(cs2)#

    #Create polyhedron for flow
    p_flow = ppl.C_Polyhedron(len(new_var_list),'universe')
    if flow_constraint_list != True:
        flow_cs = create_constraints_flow(flow_constraint_list,var_list,new_var_list)
        p_flow.add_constraints(flow_cs)
        inv_cs = create_constraints_inv(inv_constraint_list, var_list, new_var_list)#
        p_flow.add_constraints(inv_cs)#

    #combine the polyhedra
    final_poly = p1
    final_poly.concatenate_assign(p2)
    final_poly.add_space_dimensions_and_embed(1)
    final_poly.intersection_assign(p_flow)
    #print(final_poly.constraints())

    #add additional constraints
    final_poly.add_constraints(additional_cs)

    if final_poly.is_empty(): return 0
    target = final_poly.maximize(eval(expr))
    if target['bounded']: return 1
    else: return -1

#A basic linear program for calculating weight of an edge between two polyhedra
def LP(i,j,alpha,beta,cs1,cs2,flow_constraint_list,inv_constraint_list,var_list):
    new_var_list = []
    for var in var_list:
        new_var_list.append(var+"1")
    for var in var_list:
        new_var_list.append(var+"2")
    new_var_list.append("t")
    #print(new_var_list)
    declare_vars(new_var_list)

    #Create polyhedron for cs1
    p1 = ppl.C_Polyhedron(len(var_list),'universe')#
    p1.add_constraints(cs1)

    #Create polyhedron for cs2
    p2 = ppl.C_Polyhedron(len(var_list),'universe')#
    p2.add_constraints(cs2)#

    #Create polyhedron for flow
    p_flow = ppl.C_Polyhedron(len(new_var_list),'universe')
    if flow_constraint_list != True:
        flow_cs = create_constraints_flow(flow_constraint_list,var_list,new_var_list)
        p_flow.add_constraints(flow_cs)
        inv_cs = create_constraints_inv(inv_constraint_list, var_list, new_var_list)#
        p_flow.add_constraints(inv_cs)#

    #combine the polyhedra
    final_poly = p1
    final_poly.concatenate_assign(p2)
    final_poly.add_space_dimensions_and_embed(1)
    final_poly.intersection_assign(p_flow)
    #print(final_poly.constraints())

    #add additional constraints
    additional_cs = ppl.Constraint_System()
    additional_cs.insert(eval("t>=0"))
    for k in range(len(var_list)):
        if k != j:
            ineq1 = new_var_list[k] + "<=1"
            additional_cs.insert(eval(ineq1))
            ineq2 = new_var_list[k] + ">=-1"
            additional_cs.insert(eval(ineq2))
        else:
            ineq = new_var_list[k] + "==" + str(beta)
            additional_cs.insert(eval(ineq))
    ineq = str(alpha) + "*" + new_var_list[i+len(var_list)] + ">=0"
    additional_cs.insert(eval(ineq))
    final_poly.add_constraints(additional_cs)
    #print(final_poly.is_empty())
    #print(final_poly.constraints())

    expr_to_max = str(alpha) + "*" + new_var_list[i+len(var_list)]
    if check_explosion(cs1,additional_cs,expr_to_max,flow_constraint_list,\
        inv_constraint_list,var_list) == 0: 
        return 0
    elif check_explosion(cs1,additional_cs,expr_to_max,flow_constraint_list,\
        inv_constraint_list,var_list) == -1:
        return None
    else:
        if final_poly.is_empty(): return 0
        else:
            #print("%d %d %d %d "%(i,j,alpha,beta),end="")
            #print(target)
            target = final_poly.maximize(eval(expr_to_max))
            if target['bounded']: return gmpy2.mpq(target['sup_n'],target['sup_d'])
            else: return None

'''
x0 = ppl.Variable(0)
x1 = ppl.Variable(1)
x2 = ppl.Variable(2)
cs1 = ppl.Constraint_System()
cs1.insert(x0==0)
cs1.insert(x1>=0)
cs1.insert(x0+x1>=0)
cs2 = ppl.Constraint_System()
cs2.insert(x1==0)
cs2.insert(x0>=0)
cs2.insert(x0+x1>=0)
flow_list = [{'y':1,'x':1,'const':3,'sgn':'='}]
print(LP(0,1,-1,1,cs1,cs2,flow_list,['x','y','z']))
#>>0
'''

#Find weight of edge between two polyhedra
def calculate_edge_weight_PPL(cs1,cs2,flow_constraint_list,inv_constraint_list,var_list):
    value_list = []
    for i in range(len(var_list)):
        for j in range(len(var_list)):
            for alpha in {-1,1}:
                for beta in {-1,1}:
                    value_list.append(LP(i,j,alpha,beta,cs1,cs2,\
                        flow_constraint_list,inv_constraint_list,var_list))
                    #print("%d %d %d %d "%(i,j,alpha,beta),end="")
                    #print(value_list[-1])

    #print(value_list)
    if None in value_list: return None
    else: return max(value_list)

'''
x0 = ppl.Variable(0)
x1 = ppl.Variable(1)
cs1 = ppl.Constraint_System()
cs1.insert(x0>=0)
cs1.insert(x1==0)
cs2 = ppl.Constraint_System()
cs2.insert(x1>=0)
cs2.insert(x0==0)
flow_list = [{'x1':1,'x0':1,'const':0,'sgn':'='}]
inv_list = [{'x0':1,'const':0,'sgn':'>='},{'x1':1,'const':0,'sgn':'>='}]
print(calculate_edge_weight_PPL(cs1,cs2,flow_list,inv_list,['x0','x1']))
#>>1
'''

#create a dictionary
#for each facet store the two partitions it belongs to
def create_facet_part_corres(facet_list):
    corres_dict = {}

    for cs in facet_list:
        part_list = []
        cs1 = ppl.Constraint_System()
        cs2 = ppl.Constraint_System()
        for constraint in cs:
            if constraint.is_equality():
                coeff_l = list(constraint.coefficients())
                const = constraint.inhomogeneous_term()
                ineq1 = ppl.Linear_Expression(coeff_l,const)>=0
                ineq2 = ppl.Linear_Expression(coeff_l,const)<=0
                cs1.insert(ineq1)
                cs2.insert(ineq2)
            else:
                cs1.insert(constraint)
                cs2.insert(constraint)
        part_list.append(ppl.C_Polyhedron(cs1))
        part_list.append(ppl.C_Polyhedron(cs2))
        corres_dict[str(cs)] = part_list

    return corres_dict

'''
expr_list = readExprFile("equations1.txt")
#print(expr_list)
l = createAllFacetsPPL(expr_list,['x','y'])
facet_dict = create_facet_part_corres(l)
for facet_str in facet_dict:
    print("%s : "%facet_str,end="")
    print(facet_dict[facet_str][0].constraints(),end="")
    print(", ",end="")
    print(facet_dict[facet_str][1].constraints())
'''
#>>Constraint_System {x0==0, -x1>=0, -x0-x1>=0} : Constraint_System {x0>=0, -x1>=0, -x0-x1>=0}, Constraint_System {-x0>=0, -x1>=0, -x0-x1>=0}
#Constraint_System {x0==0, x1>=0, x0+x1>=0} : Constraint_System {x0>=0, x1>=0, x0+x1>=0}, Constraint_System {-x0>=0, x1>=0, x0+x1>=0}
#Constraint_System {x1==0, -x0>=0, -x0-x1>=0} : Constraint_System {x1>=0, -x0>=0, -x0-x1>=0}, Constraint_System {-x1>=0, -x0>=0, -x0-x1>=0}
#Constraint_System {x1==0, x0>=0, x0+x1>=0} : Constraint_System {x1>=0, x0>=0, x0+x1>=0}, Constraint_System {-x1>=0, x0>=0, x0+x1>=0}
#Constraint_System {x0+x1==0, -x0>=0, x1>=0} : Constraint_System {x0+x1>=0, -x0>=0, x1>=0}, Constraint_System {-x0-x1>=0, -x0>=0, x1>=0}
#Constraint_System {x0+x1==0, x0>=0, -x1>=0} : Constraint_System {x0+x1>=0, x0>=0, -x1>=0}, Constraint_System {-x0-x1>=0, x0>=0, -x1>=0}

#check if there is non null intersection between 2 lists of polyhedra
def check_intersect_polylist(poly_l1,poly_l2):
    for poly in poly_l1:
        if poly in poly_l2: return True
    return False

#create the abstraction graph
def create_abs_graph(PPCD):
    declare_vars(PPCD["Vars"])

    G = nx.DiGraph()
    facet_list = []

    node_dict_main = {}
    node_dict_aux = {}
    for edge in PPCD["Edges"]:
        temp = PPCD["Edges"][edge]
        facet_cs = ppl.Constraint_System()
        for eq_dict in temp[2]:
            eq_sgn = eq_dict.pop("sgn")
            if eq_sgn == "=":
                ineq = build_ineq_from_expr(eq_dict, "==")
            else:
                ineq = build_ineq_from_expr(eq_dict, eq_sgn)
            facet_cs.insert(eval(ineq))
        node_dict_main[temp[0]] = facet_cs
        for key in temp[1]:
            node_dict_aux[key] = facet_cs
        facet_list.append(facet_cs)

    for edge in PPCD["Edges"]:
        #print(PPCD["Edges"][edge])
        temp = PPCD["Edges"][edge]
        node1 = (temp[0],str(node_dict_main[temp[0]]))
        for key in temp[1]:
            node2_1 = (key,str(node_dict_main[key]))
            node2_2 = (key,str(node_dict_aux[key]))
            weight = calculate_edge_weight_PPL(node_dict_main[temp[0]],node_dict_main[key],\
                PPCD["Flow"][key],PPCD["Inv"][key],PPCD["Vars"])
            
            if weight == None:
                G.add_node(node1,type = "r")
                G.add_node(node2_2,type = "r")
                G.add_edge(node1,node2_2,type="p",prob=temp[1][key],weight=weight)
            else:
                G.add_node(node1,type = "r")
                G.add_node(node2_1,type = "r")
                G.add_edge(node1,node2_1,type="p",prob=temp[1][key],weight=weight)
    
    #for node1,node2,data in G.edges.data(): print(data["weight"])
    return G,facet_list
    
#draw the abstraction graph
def draw_abs_graph(G,facet_list):
    facet_dict = {}
    for facet in facet_list:
        facet_dict[str(facet)] = "f" + str(facet_list.index(facet))
    #print(facet_dict)

    color_map = []
    node_size = []
    edge_labels = {}
    node_labels = {}

    for node,data in G.nodes.data():
        if data['type'] == 'r':
            color_map.append('pink')
            node_size.append(1000)
            node_labels[node] = node[0] + "," + facet_dict[node[1]]
        else:
            color_map.append('green')
            node_size.append(20)

    for node1,node2,data in G.edges.data():
        if data['type'] == 'c':
            edge_labels[(node1,node2)] = data['weight']
        elif data['type'] == 'p':
            edge_labels[(node1,node2)] = str(data['prob']) + "," +\
                str(data['weight'])

    pos = nx.spring_layout(G)
    plt.figure()    
    nx.draw(G,pos,node_size=node_size,node_color=color_map,\
        labels=node_labels,font_size=5)
    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels,\
        font_color='black')
    plt.axis('off')
    plt.show()

'''
PPCD = ppcd_input.readPPCD("ppcd_example3.txt")
G,facet_list = create_abs_graph(PPCD)
#print(list(G.edges().data()))
draw_abs_graph(G,facet_list)
#print(list(G.nodes()))
'''

#Convert graph G to a graph to be analyzed
#source = (loc,str(facet))
#keep only edges reachable from source
#weight = 0 if edge type p or d
#weight = -log(weight)/None if egde type c
#if a cycle had weight>1 in the new graph it has weight<0
def graph_for_analysis(G,source):
    G_anal = nx.DiGraph()
    reachable_edges = list(nx.edge_dfs(G,source,orientation='original'))
    for e in reachable_edges:
        #print("%s,%s"%(e[0],e[1]))
        edge_data = G.get_edge_data(e[0],e[1])
        if edge_data['type'] == 'c':
            G_anal.add_node(e[0],type='r')
            G_anal.add_node(e[1],type='r')
            if edge_data['weight'] == None:
                G_anal.add_edge(e[0],e[1],weight=None)
            else:
                G_anal.add_edge(e[0],e[1],weight=-(math.log(edge_data['weight'])))
        elif edge_data['type'] == 'p':
            G_anal.add_node(e[0],type='r')
            G_anal.add_node(e[1],type='r')
            if edge_data['weight'] == None:
                G_anal.add_edge(e[0],e[1],weight=None,prob=edge_data['prob'])
            elif edge_data['weight'] == 0:#
                G_anal.add_edge(e[0],e[1],weight=-math.inf,prob=edge_data['prob'])#
            else:
                G_anal.add_edge(e[0],e[1],weight=-(math.log(edge_data['weight'])),\
                    prob=edge_data['prob'])
        else:
            G_anal.add_node(e[0],type='r')
            G_anal.add_node(e[1],type='d')
            G_anal.add_edge(e[0],e[1],weight=0)

    return G_anal

'''
PPCD = ppcd_input.readPPCD("ppcd_example1.txt")
G,facet_list = create_abs_graph(PPCD)
source = list(G.nodes())[0]
G_anal = graph_for_analysis(G, source)
'''

def draw_Ganal(G,facet_list):
    facet_dict = {}
    for facet in facet_list:
        facet_dict[str(facet)] = "f" + str(facet_list.index(facet))
    #print(facet_dict)

    color_map = []
    node_size = []
    edge_labels = {}
    node_labels = {}

    for node,data in G.nodes.data():
        if data['type'] == 'r':
            color_map.append('pink')
            node_size.append(1000)
            node_labels[node] = node[0] + "," + facet_dict[node[1]]
        else:
            color_map.append('green')
            node_size.append(20)

    for node1,node2,data in G.edges.data():
        edge_labels[(node1,node2)] = data['prob'],data['weight']

    pos = nx.circular_layout(G)
    plt.figure()    
    nx.draw(G,pos,node_size=node_size,node_color=color_map,\
        labels=node_labels,font_size=5)
    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels,\
        font_color='black')
    plt.axis('off')
    plt.show()


'''
PPCD = ppcd_input.readPPCD("ppcd_example.txt")
G,facet_list = create_abs_graph(PPCD)
source = list(G.nodes())[0]
G_anal = graph_for_analysis(G, source)
draw_Ganal(G_anal, facet_list)
'''

