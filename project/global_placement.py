"""# Simulated Annealing

"""

import math
import random
import numpy as np
import LEFDEFParser
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider

def accept(delC, T):
  if delC <= 0: return True
  return random.random() < math.exp(-delC/T)

# costEval(F): cost evaluation routine
# accept(\Delta C, T): accept routine
# perturb(F): perturbation routine
def sa(Tmin, Tmax, N, alpha, V, Vdict, Pdict, bb, init_cost, plot = True ):
  assert(alpha < 1. and Tmin < Tmax)
  T = Tmax
  C = init_cost

  minC = C
  minV = V[:]
  Clist = []
  Temp = []
  wT = bb[1][0] - bb[0][0]
  hT = bb[1][1] - bb[0][1]

  while T > Tmin:
    for i in range(N):
      savedSol = [np.array([v._comp.location().x for v in V]), np.array([v._comp.location().y for v in V])]
      Vnew, v, (v1,v2), p_case = perturb(V, wT, hT, bb)
      Cnew = costUpdate(v, (v1,v2), p_case, Vnew, Vdict, Pdict, savedSol, C)

      if accept(Cnew - C, T):
        C, V = Cnew, Vnew
        print("Current cost = ", C)
        if minC >= Cnew:
          minC, minV = Cnew, Vnew
        Clist.append(Cnew)
        Temp.append(T)
      else:
        for i in range(len(V)):
          V[i]._comp.setLocation(savedSol[0][i], savedSol[1][i])

    wT = wT * (math.log(T * alpha) / math.log(T))
    hT = hT * (math.log(T * alpha) / math.log(T))
    T = T * alpha

  if plot:
    import matplotlib.pyplot as plt
    plt.plot(Temp, Clist)
    plt.xlim(max(Temp), min(Temp))
    plt.xscale('log')
  return minV, minC

"""# Perturb function


"""

# Class to hold neighbours and pins connected to each standard cell
class Vertex:
    def __init__(self, c, dim):
        self._comp = c
        self._index = -1
        self._dim = dim
        self._nets = []

    def __repr__(self):
        return self._comp.name() + f' {self._nbrs} {self._pinNbrs}'

def perturb(Vi, wT, hT, bb):
    done = False
    p_case = -1
    bb_ll, bb_ur = bb[0], bb[1]
    while(done == False):
      v1 = -1
      v2 = -1
      v = -1
      choice = random.randint(1,2)
      match choice:
        case 1: # MOVE: Shift a cell to a new position
          v = random.choice(Vi)
          pos = (v._comp.location().x + v._dim[0]/2, v._comp.location().y + v._dim[1]/2)
          x_new_ll =  pos[0] + random.uniform(-wT/2, wT/2) - v._dim[0]/2
          y_new_ll =  pos[1] + random.uniform(-hT/2, hT/2) - v._dim[1]/2
          if (x_new_ll > bb_ll[0] and y_new_ll > bb_ll[1] and x_new_ll + v._dim[0] < bb_ur[0] and y_new_ll + v._dim[1] < bb_ur[1]):
            v._comp.setLocation(int(x_new_ll), int(y_new_ll))
            done = True
            p_case = 1

        case 2: # SWAP: Exchange two cells
          v1 = random.choice(Vi)
          pos1 = (v1._comp.location().x + v1._dim[0]//2, v1._comp.location().y + v1._dim[1])

          i = 0
          iter_limit = 20
          while(True):
            i += 1
            if(i>= iter_limit):
              break
            v2 = random.choice(Vi)
            pos2 = (v2._comp.location().x + v2._dim[0]//2, v2._comp.location().y + v2._dim[1])

            x_diff = abs(pos1[0] - pos2[0])
            y_diff = abs(pos1[1] - pos2[1])

            if(x_diff <= wT and y_diff <= hT and v1 != v2):
              v1_new_ll = (v2._comp.location().x, v2._comp.location().y)
              v2_new_ll = (v1._comp.location().x, v1._comp.location().y)

              if(v1_new_ll[0] > bb_ll[0] and v1_new_ll[1] > bb_ll[1] and v1_new_ll[0] + v1._dim[0] < bb_ur[0] and v1_new_ll[1] + v1._dim[1] < bb_ur[1]):
                if(v2_new_ll[0] > bb_ll[0] and v2_new_ll[1] > bb_ll[1] and v2_new_ll[0] + v2._dim[0] < bb_ur[0] and v2_new_ll[1] + v2._dim[1] < bb_ur[1]):
                  temp = (v1._comp.location().x, v1._comp.location().y)
                  v1._comp.setLocation(v2._comp.location().x, v2._comp.location().y)
                  v2._comp.setLocation(temp[0], temp[1])
                  break

          done = True
          p_case = 2

    return Vi, v, (v1,v2), p_case

"""# Cost Eval"""

# Calculate overlap between two cells
def overlap(p1, p2):
    xd = max(min(p1[2] - p2[0], p2[2] - p1[0]), 0)
    yd = max(min(p1[3] - p2[1], p2[3] - p1[1]), 0)
    return xd * yd * 1.e-6

def overlapCost(V):
    overlapArea = 0.
    P = [(v._comp.location().x, v._comp.location().y, v._comp.location().x + v._dim[0], v._comp.location().y + v._dim[1]) for v in V]
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            overlapArea += overlap(P[i], P[j])
    print("overlap area = ", overlapArea)

    return overlapArea

def hpwlCost(nets, Vdict, Pdict, V):
    hpwl = 0.
    for n in nets:
        minx, miny = 1.e20, 1.e20
        maxx, maxy = -1.e20, -1.e20
        for p in n.pins():
            if p[0] != 'PIN':
                assert(p[0] in Vdict)
                v = V[Vdict[p[0]]]
                pos = (v._comp.location().x + v._dim[0]//2, v._comp.location().y + v._dim[1])
            else:
                assert(p[1] in Pdict)
                pos = (Pdict[p[1]].x, Pdict[p[1]].y)
            minx = min(minx, pos[0])
            maxx = max(maxx, pos[0])
            miny = min(miny, pos[1])
            maxy = max(maxy, pos[1])
        hpwl += ((maxx - minx) + (maxy - miny)) * 1.e-3
    return hpwl

def costEval(V, Vdict, Pdict, d):
  oCost = overlapCost(V)
  hCost = hpwlCost(d.nets(), Vdict, Pdict, V)
  total_cost = oCost + hCost
  print("overlap cost :", oCost, " hpwlCost :", hCost, "total cost :", oCost + hCost)
  return total_cost

def overlapUpdate(v, V, savedSol):
  P_old = (savedSol[0][v._index], savedSol[1][v._index], savedSol[0][v._index] + v._dim[0], savedSol[1][v._index] + v._dim[1])
  P_new = (v._comp.location().x, v._comp.location().y, v._comp.location().x + v._dim[0], v._comp.location().y + v._dim[1])
  P = [(u._comp.location().x, u._comp.location().y, u._comp.location().x + u._dim[0], u._comp.location().y + u._dim[1]) for u in V]

  overlapArea_old, overlapArea_new = 0, 0

  for i in range(len(V) - 1):
    if(i != v._index):
      overlapArea_old += overlap(P_old, P[i])
      overlapArea_new += overlap(P_new, P[i])

  return overlapArea_new - overlapArea_old

def hpwlUpdate(v, V, Vdict, Pdict, savedSol):
  hpwl_old, hpwl_new = 0, 0
  nets = v._nets

  for n in nets:
    minx, miny = 1.e20, 1.e20
    maxx, maxy = -1.e20, -1.e20
    minx_new, miny_new = 1.e20, 1.e20
    maxx_new, maxy_new = -1.e20, -1.e20
    for p in n.pins():
        if p[0] != 'PIN':
            assert(p[0] in Vdict)
            u = V[Vdict[p[0]]]
            pos = (u._comp.location().x + u._dim[0]//2, u._comp.location().y + u._dim[1]//2)
            pos_new = (u._comp.location().x + u._dim[0]//2, u._comp.location().y + u._dim[1]//2)

            if (u == v):
              pos =  (savedSol[0][v._index] + u._dim[0]//2, savedSol[1][v._index] + u._dim[1]//2)
        else:
            pos = (Pdict[p[1]].x, Pdict[p[1]].y)
            pos_new = pos
        minx = min(minx, pos[0])
        maxx = max(maxx, pos[0])
        miny = min(miny, pos[1])
        maxy = max(maxy, pos[1])
        minx_new = min(minx_new, pos_new[0])
        maxx_new = max(maxx_new, pos_new[0])
        miny_new = min(miny_new, pos_new[1])
        maxy_new = max(maxy_new, pos_new[1])
    hpwl_old += ((maxx - minx) + (maxy - miny)) * 1.e-3
    hpwl_new += ((maxx_new - minx_new) + (maxy_new - miny_new)) * 1.e-3

  return hpwl_new - hpwl_old

def costUpdate(v, v_swap, p_case, V, Vdict, Pdict, savedSol, old_cost):
  match p_case:
    case 1: # MOVE: Shift a cell to a new position
      o_delta = overlapUpdate(v, V, savedSol)
      h_delta = hpwlUpdate(v, V, Vdict, Pdict, savedSol)
    case 2: # SWAP: Exchange two cells
      o_delta_1 = overlapUpdate(v_swap[0], V, savedSol)
      o_delta_2 = overlapUpdate(v_swap[1], V, savedSol)
 
      o_delta = o_delta_1 + o_delta_2

      h_delta_1 = hpwlUpdate(v_swap[0], V, Vdict, Pdict, savedSol)
      h_delta_2 = hpwlUpdate(v_swap[1], V, Vdict, Pdict, savedSol)
      h_delta = h_delta_1 + h_delta_2

    case _:
      print("Wrong case from perturb")

  new_cost = old_cost + o_delta + h_delta

  return new_cost

"""# Plotter"""

def plot(sol, bb):

    fig, ax = plt.subplots()
    ax.plot(sol[0][0], sol[0][1], 'o')
    ax.set_xlim([bb[0][0], bb[1][0]])
    ax.set_ylim([bb[0][1], bb[1][1]])
    fig.subplots_adjust(bottom=0.25)

    plt.show()

# load the DEF file and build the connectivity graph using the Vertex class
# Boundary pins have the name PIN followed by the pinName
def place(deffile, leffile, outfile):
    # load LEF file for cell dimensions
    l = LEFDEFParser.LEFReader()
    dimLookup = dict()
    if leffile:
      l.readLEF(leffile)
      dimLookup = {m.name():(m.xdim(), m.ydim()) for m in l.macros()}

    d = LEFDEFParser.DEFReader()
    d.readDEF(deffile)
    chip_bbox = d.bbox()
    V = [Vertex(c, dimLookup[c.macro()]) for c in d.components()]
    for i in range(len(V)):
        V[i]._index = i
    Vdict = {V[i]._comp.name():i for i in range(len(V))}

    P = [p.origin() for p in d.pins()]
    Pdict = {d.pins()[i].name():P[i] for i in range(len(P))}

    for n in d.nets():
        u = list()
        pins = list()
        for p in n.pins():
            if p[0] != 'PIN':
                assert(p[0] in Vdict)
                u.append(Vdict[p[0]])
            else:
                assert(p[1] in Pdict)
                pins.append(Pdict[p[1]])
        for i in range(len(u)):
            V[u[i]]._nets.append(n)

    bb = ((chip_bbox.ll.x, chip_bbox.ll.y),(chip_bbox.ur.x, chip_bbox.ur.y))

    for v in V:
      v._comp.setLocation(random.randint(bb[0][0], bb[1][0] - v._dim[0]), random.randint(bb[0][1], bb[1][1] - v._dim[1]))

    initial_coords = [(np.array([v._comp.location().x for v in V]), np.array([v._comp.location().y for v in V]))]
    plot(initial_coords, bb)

    Tmax = 1000
    oCost = overlapCost(V)
    hCost = hpwlCost(d.nets(), Vdict, Pdict, V)
    total_cost = oCost + hCost
    print("overlap cost :", oCost, " hpwlCost :", hCost, "total cost :", oCost + hCost)

    Vmin, Cmin = sa(10, Tmax, 10*len(V), 0.9, V, Vdict, Pdict, bb, total_cost, True)

    final_coords = [(np.array([v._comp.location().x for v in Vmin]), np.array([v._comp.location().y for v in Vmin]))]

    plot(final_coords, bb)
    d.writeDEF(outfile)
    oCost = overlapCost(V)
    hCost = hpwlCost(d.nets(), Vdict, Pdict, V)
    print("overlap cost :", oCost, " hpwlCost :", hCost, "total cost :", oCost + hCost)

place('dma.def', 'dma.lef', 'dma_out.def')