import random
import gmpy2

def CreateDist(loclist):
    l = [10*random.randint(1,10) for i in range(len(loclist))]
    s = sum(l)
    dist = [gmpy2.mpq(element,s) for element in l]
    dist_dict = {}
    for i in range(len(loclist)):
        dist_dict[loclist[i]] = dist[i]

    return dist_dict
  
#print(CreateDist(['q1','q2','q3','q4','q5']))

def CreateLocs2D(div_num,state_num):
    Q = []
    for i in range(4):
        Q.append([])
        for j in range(div_num):
            Q[i].append([])
            for k in range(state_num):
                Q[i][j].append('q'+str(i)+str(j)+str(k))
    
    return Q

#print(CreateLocs2D(2,2))

def CreateFlow2D(Q,div_num,state_num,Vars):
    Flow = {}
    for i in range(4):
        for j in range(div_num):
            for k in range(state_num):
                rate = random.randint(1,5)
                if i==0: Flow[Q[i][j][k]] = \
                    Vars[0]+"+"+str(rate)+Vars[1]+"=0"
                elif i==1: Flow[Q[i][j][k]] = \
                    "-"+str(rate)+Vars[0]+"+"+Vars[1]+"=0"
                elif i==2: Flow[Q[i][j][k]] = \
                    "-"+Vars[0]+"-"+str(rate)+Vars[1]+"=0"
                else: Flow[Q[i][j][k]] = \
                    str(rate)+Vars[0]+"-"+Vars[1]+"=0"

    return Flow

#Q = CreateLocs2D(2,2)
#print(CreateFlow2D(Q,2,2,["x","y"]))

def CreateEdge2D(Q,div_num,state_num,guard_list):
    edges = {}
    num = 1
    for i in range(4):
        for j in range(div_num):
            for k in range(state_num):
                if j==div_num-1 and i==3:
                    edges['e'+str(num)] = \
                        Q[i][j][k],CreateDist(Q[0][0]),guard_list[i][j]
                elif j==div_num-1 and i<3:
                    edges['e'+str(num)] = \
                        Q[i][j][k],CreateDist(Q[i+1][0]),guard_list[i][j]
                else:
                    edges['e'+str(num)] = \
                        Q[i][j][k],CreateDist(Q[i][j+1]),guard_list[i][j]
                num += 1

    return edges

#Q = CreateLocs2D(2,2)
#print(CreateEdge2D(Q,2,2))

'''
def CreateReset(edges):
    reset = {}
    for edge in edges:
        for key in edges[edge][1]:
            reset[edge,key] = True
    
    return reset

#Q = CreateLocs2D(2,2)
#edges = CreateEdge2D(Q,2,2)
#print(CreateReset(edges))
'''

def WriteFile(fname,Q,div_num,state_num,Vars,Flow,edges,inv_list,guard_list,self_prob):
    fp = open(fname,"w")
    #print Q
    print("#Locations{",file=fp)
    loc_str = ""
    inv_str = ""
    for i in range(4):
        for j in range(div_num):
            for k in range(state_num):
                loc_str += Q[i][j][k]+","
                inv_str += Q[i][j][k]+" : "+inv_list[i][j]+"\n"
    print(loc_str[:-1],file=fp)
    print("}\n",file=fp)

    #print Vars
    print("#Dimensions{",file=fp)
    var_str = ""
    for dim in Vars: var_str += dim+","
    print(var_str[:-1],file=fp)
    print("}\n",file=fp)

    #print Inv
    print("#Inv{",file=fp)
    print(inv_str[:-1],file=fp)
    print("}\n",file=fp)

    #print Flow
    print("#Flow{",file=fp)
    for q in Flow:
        print("%s : %s"%(q,Flow[q]),file=fp)
    print("}\n",file=fp)

    #print edges
    print("#Edges{",file=fp)
    for e in edges:
        print("%s : "%e,end="",file=fp)
        print("%s,"%edges[e][0],end="",file=fp)
        print("(",end="",file=fp)
        dist_str = edges[e][0]+"="+str(self_prob)+","
        for q in edges[e][1]:
            dist_str += q+"="+str(edges[e][1][q]*(1-self_prob))+","
        print(dist_str[:-1],end="",file=fp)
        print("),",end="",file=fp)
        print(edges[e][2],file=fp)
    print("}\n",file=fp)

    '''
    #print reset
    print("#Reset{",file=fp)
    for e,q in reset:
        print("%s,%s : True"%(e,q),file=fp)
    print("}",file=fp)
    '''

    fp.close()


Q = CreateLocs2D(2,12)
Vars = ['x','y']
inv_list = [["x>=0,y>=0,x-y>=0","x>=0,y>=0,x-y<=0"],\
    ["x<=0,y>=0,x+y>=0","x<=0,y>=0,x+y<=0"],\
        ["x<=0,y<=0,x-y<=0","x<=0,y<=0,x-y>=0"],\
            ["x>=0,y<=0,x+y<=0","x>=0,y<=0,x+y>=0"]]
guard_list = [["(x-y=0,y>=0,x>=0)","(x=0,y>=0)"],\
    ["(x<=0,y>=0,x+y=0)","(x<=0,y=0)"],\
        ["(x-y=0,y<=0,x<=0)","(x=0,y<=0)"],\
            ["(x>=0,y<=0,x+y=0)","(x>=0,y=0)"]]

Flow = CreateFlow2D(Q,2,12,Vars)
edges = CreateEdge2D(Q,2,12,guard_list)
#reset = CreateReset(edges)
WriteFile("ppcd_example4.txt",Q,2,12,Vars,Flow,edges,inv_list,guard_list,gmpy2.mpq(5,100))

