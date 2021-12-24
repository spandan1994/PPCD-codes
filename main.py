import ppcd_input
import ppcd_abstract
import check_stability
import timeit
from statistics import mean

'''
def print_square(x):
    return (x**2)

t = timeit.repeat(lambda: print_square(3), number=10, repeat=5)
  
# printing the execution time
print(t)
'''

filename = input("Provide PPCD file name : ")
start = timeit.default_timer()
PPCD = ppcd_input.readPPCD(filename)
G,facet_list = ppcd_abstract.create_abs_graph(PPCD)
end = timeit.default_timer()
print(list(G.nodes()))
nodenum = eval(input("Provide source node number : "))
source = list(G.nodes())[nodenum]
#print(G.edges.data())
#ppcd_abstract.draw_abs_graph(G,facet_list)

print("Time to create abstract graph : %f"%(end-start))

as_t_list = timeit.repeat(lambda: check_stability.check_abs_stability(G, source), \
    number=1, repeat=5)
r = check_stability.check_abs_stability(G, source)
if r==-1:
    print("Exploding region reachable from source")
elif r==1:
    print("System is absolutely stable")
else:
    print("System is not absolutely stable")
print("Time of completion : %f"%mean(as_t_list))

as_t_list = timeit.repeat(lambda: check_stability.check_as_stability(G, source), \
    number=1, repeat=5)
r = check_stability.check_as_stability(G, source)
if r==1:
    print("System is almost surely stable")
else:
    print("System is not almost surely stable")
print("Time of completion : %f"%mean(as_t_list))
