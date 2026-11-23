# Name: Srichandra Siripuram
# Roll No: EE20B128
# Course number: EE5333
# Assignment: Graph colouring and bipartioning using ILP


# N is num. vertices; vertices are indexed {0,1,...,(N-1)}
# E is a list of edges; each edge is an unordered pair of vertices
# C is list of color costs : C[i] is the cost of color i
#                            color ranges from [0,1,...,len(C)-1]
# Return value : list of colors if colorable empty list otherwise;
#                ith color in the list is for the ith vertex

##### CODE FOR GRAPH COLOURING #####

def colorGraph(N, E, C):
  import mip
  model = mip.Model("GraphColouring")
  model.verbose = 0 # turn off Cbc console logs
  x = []
  for v in range(N):
    x.append([model.add_var(var_type=mip.BINARY,\
        name=f"x_{v}_{i}") for i in range(len(C))])
  
  model.objective = mip.minimize(mip.xsum(C[c]*x[i][c] for i in range(N) for c in range(len(C))))

  # constraint 1
  for i in range(N):
    model += mip.xsum(x[i]) == 1

  # constraint 2
  for e in E:
    for i in range(len(C)):
      model += (x[e[0]][i] + x[e[1]][i] <= 1)

  model.optimize()
  model.write("graph_colouring.lp")
  if model.status == mip.OptimizationStatus.OPTIMAL:
    return [c for v in range(N) for c in range(len(C)) if x[v][c].x >= 0.9]
  return []


##### CODE FOR BIPARTITION-2 #####

def bipartition2(N, E):
  import mip
  model = mip.Model("Bi-partition2")

  # variables for vertex 'v' in partition A or partition B
  x_v_a = [model.add_var(f"x_{v}_a", var_type = mip.BINARY) for v in range(N)]
  x_v_b = [model.add_var(f"x_{v}_b", var_type = mip.BINARY) for v in range(N)]

  # variables for an edge in partition A or partition B
  x_uv_a = [model.add_var(f"x_{u}_{v}_a", var_type = mip.BINARY) for u,v in E]
  x_uv_b = [model.add_var(f"x_{u}_{v}_b", var_type = mip.BINARY) for u,v in E]
  model.verbose = 0

  model.objective = mip.maximize(mip.xsum(x_uv_a + x_uv_b))

  # constraint 1
  model += (mip.xsum(x_v_a) == N//2)
  model += (mip.xsum(x_v_b) == N//2)

  # constraint 3, 4
  for e, (u,v) in enumerate(E):
    model += x_uv_a[e] <= x_v_a[u]
    model += x_uv_a[e] <= x_v_a[v]
    model += x_uv_b[e] <= x_v_b[u]
    model += x_uv_b[e] <= x_v_b[v]

  # constraint 5
  for i in range(N):
    model += x_v_a[i] + x_v_b[i] == 1 

  model.write("bipartition2.lp")
  model.optimize()
  if model.status == mip.OptimizationStatus.OPTIMAL:
    A = [i for i in range(N) if x_v_a[i].x >= 0.9]
    B = [i for i in range(N) if x_v_b[i].x >= 0.9]
    return (A, B, len(E) - model.objective.x) # returning the partitions A, B and the number of cuts
  return None