import re
import networkx as nx
import gmpy2
import matplotlib.pyplot as plt

'''
Read a linear equation and return as dictionary
Example : x+y<=3 will be returned as 
{'const': 3.0, 'sgn': '<=', 'x': 1.0, 'y': 1.0}
'''
def readEq(equation):
    d = {"const":0.0,"sgn":"="}
    tempeqn = equation
    pattern1 = "\s*([0-9+-.]*)([a-zA-Z]+[0-9_']*)\s*" #identify var
    pattern2 = "\s*([><=]+)\s*([0-9+-.]+)\s*" #identify const with sgn
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
            tempeqn = re.sub(pattern1,"",tempeqn,1) #replace identified const with null str
        except:
            result = re.search(pattern2,tempeqn)
            d["const"] = (eval(result.group(2)))*1.0
            d["sgn"] = result.group(1)
            break
    
    return d

#s = " x +y <=  3"
#print(readEq(s))

'''
Reads a probability distribution and returns a dictionary
Example : 
q1=0.3, q2=0.7 as
{'q1': 0.3, 'q2': 0.7}
'''
def readDist(distribution):
    d = {}
    tempstr = distribution
    pattern = "\s*(\w+)\s*[=]\s*(([0-9.]+)[/]([0-9]+))\s*[,]*" #identify one key-val
    while tempstr != "":
        result = re.search(pattern,tempstr)
        value = gmpy2.mpq(eval(result.group(3)),eval(result.group(4))) #evaluate probability
        d[result.group(1)] = value
        tempstr = re.sub(pattern,"",tempstr,1) #remove key-val

    return d

#s = "q1= 3/70, q2 =67/70 "
#print(readDist(s))

'''
Reads locations given a string and outputs list
Example : q1,q2,q3 as ['q1','q2','q3']
'''
def readLocs(locations):
    l = []
    tempstr = locations
    pattern = "\s*(\w+)\s*[,]*" #identify a single location
    while tempstr != "":
        result = re.search(pattern,tempstr)
        l.append(result.group(1)) #append location
        tempstr = re.sub(pattern,"",tempstr,1) #remove identified location from str
    
    return l

#s = "q1,q2 , q3, q4"
#print(readLocs(s))

'''
Reads variables given a string and outputs list
Example : x,y,z as ['x','y','z']
'''
def readVars(variables):
    l = []
    tempstr = variables
    pattern = "\s*([a-zA-Z]+[0-9_]*)\s*[,]*" #identify a single var
    while tempstr != "":
        result = re.search(pattern,tempstr)
        l.append(result.group(1)) #append var
        tempstr = re.sub(pattern,"",tempstr,1) #remove identified var from str
    
    return l

#s = "x,y , z_1, z_2"
#print(readVars(s))

'''
Reads Inv given a string and adds output to a dictionary
Example : 
q1 : x>=0,y>=0 added as 
{'q1':[{'x':1,'sgn':'>=','const':0},{'y':1,'sgn':'>=','const':0]}
'''
def readInv(invline,d): #d is a dict
    #d = {}
    tempstr = invline
    pattern = "\s*(\w+)\s*[:]\s*(.+)" #identify Inv pattern
    result = re.search(pattern,tempstr)
    list_lin_const = []
    str_lin_const = result.group(2)
    if re.match("^True",str_lin_const): #check no constraint case
        d[result.group(1)] = True
        return d
    
    lin_pattern = "\s*([^ ,\n]+)\s*[,]*" #identify single constraint pattern
    while str_lin_const != "":
        lin_const = re.search(lin_pattern,str_lin_const)
        new_dict = readEq(lin_const.group(1)) #obtain a lin constraint
        list_lin_const.append(new_dict) #add lin constraint to list
        str_lin_const = re.sub(lin_pattern,"",str_lin_const,1) #remove lin constraint from str
    
    d[result.group(1)] = list_lin_const
    return d

'''
d = {}
s = "q1 : x>=0, y>=0 "
print(readInv(s,d))
s = "q2 : x>=0, y>=0 "
print(readInv(s,d))
s = "q3 : True"
print(readInv(s,d))
'''

#Use readInv for reading Flows as well since they have same structure
#s = "q1 : x=0.7,y=-0.7,x-z=0"
#d = {}
#print(readInv(s,d))

'''
Reads an edge and adds output to a dictionary
Example :
e1 : q1,(q2=0.3,q3=0.7),(x=0,y+0.714286x>=0) added as
{'e1': ('q1',{distribution},[guard_lin_consts])}
'''
def readEdge(uniedge,d): #d is a dict
    tempstr = uniedge
    pattern = "\s*(\w+)\s*[:]\s*(\w+)\s*[,]\s*[(]([^\n)(]+)[)]\s*[,]\s*[(]([^\n)(]+)[)]" #identify edge pattern
    result = re.search(pattern,tempstr)

    #get guard lin constraints
    list_lin_const = []
    str_lin_const = result.group(4)
    if re.match("^True",str_lin_const): #check no constraint case
        list_lin_const = True
    else:
        lin_pattern = "\s*([^ ,\n]+)\s*[,]*" #identify single constraint pattern
        while str_lin_const != "":
            lin_const = re.search(lin_pattern,str_lin_const)
            new_dict = readEq(lin_const.group(1)) #obtain a lin constraint
            list_lin_const.append(new_dict) #add lin constraint to list
            str_lin_const = re.sub(lin_pattern,"",str_lin_const,1) #remove lin constraint
        #print(list_lin_const)

    #get distribution
    new_dist = readDist(result.group(3))
    #print(new_dist)

    #build tuple
    edge_tuple = tuple((result.group(2),new_dist,list_lin_const))
    #print(edge_tuple)

    #set key-value pair
    d[result.group(1)] = edge_tuple
    return d

'''
s = "e1 : q1,(q2=0.3,q3=0.7),(x=0,y+0.714286x>=0)"
d = {}
print(readEdge(s,d))
s = "e2 : q2,(q1=0.7,q4=0.3),(x<=0,y+0.714286x=0)"
print(readEdge(s,d))
s = "e3 : q3,(q1=0.7,q4=0.3),(True)"
print(readEdge(s,d))
'''

'''
Reads a Reset condition and adds output to dict
Example :
e1,q2 : 10x'>=-5, x'<=0, y'>=-1, y'<=0 added as
{('e1','q2'): [reset_lin_consts]}
'''
'''
def readReset(resetline,d): #d is a dict
    tempstr = resetline
    pattern = "\s*(\w+)\s*[,]\s*(\w+)\s*[:]\s*(.+)"
    result = re.search(pattern,tempstr)

    list_lin_const = []
    str_lin_const = result.group(3)
    if re.match("^True",str_lin_const): #check no constraint case
        list_lin_const = True
    else:
        lin_pattern = "\s*([^ ,\n]+)\s*[,]*" #identify single constraint pattern
        while str_lin_const != "":
            lin_const = re.search(lin_pattern,str_lin_const)
            new_dict = readEq(lin_const.group(1)) #obtain a lin constraint
            list_lin_const.append(new_dict) #add lin constraint to list
            str_lin_const = re.sub(lin_pattern,"",str_lin_const,1) #remove lin const

    key = tuple((result.group(1),result.group(2)))
    d[key] = list_lin_const
    return d
'''

'''
d={}
s = "e1,q3 : True"
print(readReset(s,d))
s = "e1,q2 : 10x'>=-5, x'<=0, y'>=-1, y'<=0"
print(readReset(s,d))
'''

'''
Read a PPCD from the input file and return as dictionary
{'Locs': , 'Vars': , 'Inv': , 'Flow': , 'Edges': }
'''

def readPPCD(filename):
    fp = open(filename,"r")
    newline = fp.readline()
    #print(newline)
    while newline:
        if re.match("^#Locations",newline):
            loc_list = []
            newline = fp.readline()
            while newline != "}\n":
                loc_list = readLocs(newline)
                newline = fp.readline()
            #print(loc_list)
        elif re.match("^#Dimensions",newline):
            var_list = []
            newline = fp.readline()
            while newline != "}\n":
                var_list = readLocs(newline)
                newline = fp.readline()
            #print(var_list)
        elif re.match("^#Inv",newline):
            Inv_dict = {}
            newline = fp.readline()
            while newline != "}\n":
                Inv_dict = readInv(newline,Inv_dict)
                newline = fp.readline()
            #print(Inv_dict)
        elif re.match("^#Flow",newline):
            Flow_dict = {}
            newline = fp.readline()
            while newline != "}\n":
                Flow_dict = readInv(newline,Flow_dict)
                newline = fp.readline()
            #print(Flow_dict)
        elif re.match("^#Edges",newline):
            Edge_dict = {}
            newline = fp.readline()
            while newline != "}\n":
                Edge_dict = readEdge(newline,Edge_dict)
                newline = fp.readline()
            #print(Edge_dict)
        '''
        elif re.match("^#Reset",newline):
            Reset_dict = {}
            newline = fp.readline()
            while newline != "}\n":
                Reset_dict = readReset(newline,Reset_dict)
                newline = fp.readline()
            #print(Reset_dict)
        '''
        newline = fp.readline()

    PPCD = {}
    PPCD["Locs"] = loc_list
    PPCD["Vars"] = var_list
    PPCD["Inv"] = Inv_dict
    PPCD["Flow"] = Flow_dict
    PPCD["Edges"] = Edge_dict
    #PPCD["Reset"] = Reset_dict

    return PPCD
    
    fp.close()

#PPCD = readPPCD("ppcd_example.txt")
#print(PPCD)

#Take input PPCD as dictionary and draw corresponding graph
def graphPPCD(PPCD):
    G = nx.DiGraph()
    edge_list = [] #stores graph edges
    color_map = [] #stores a color corresponding to each node according to position
    node_size = [] #stores an integer size of node according position
    edge_labels = {} #stores edge labels
    labels = {} #stores node labels
    for edge in PPCD["Edges"]:
        #print(PPCD["Edges"][edge])
        temp = PPCD["Edges"][edge]
        edge_list.append((temp[0],str(temp[1])))
        labels[temp[0]] = temp[0]
        pi = temp[1]
        for loc in pi:
            edge_list.append((str(temp[1]),loc))
            edge_labels[(str(temp[1]),loc)] = pi[loc]
    #print(edge_list)
    #print(edge_labels)
    G.add_edges_from(edge_list)
    #print(G.nodes())
    for node in G.nodes():
        if re.match("^{",node): color_map.append('green'); node_size.append(20) #if node is a distribution put green color and small size
        else: color_map.append('pink'); node_size.append(1000) #if node is a location put pink color and larger size

    pos = nx.spring_layout(G)
    plt.figure()    
    nx.draw(G,pos,node_size=node_size,node_color=color_map,\
        labels=labels)
    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels,\
        font_color='black')
    plt.axis('off')
    plt.show()

#graphPPCD(PPCD)





