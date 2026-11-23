from queue import Queue
import copy 
import random
import math

def accept(delC, T):
  if delC <= 0: return True
  return random.random() < math.exp(-delC/T)

# F = initial slicing floorplan in Polish form
# costEval(F): cost evaluation routine
# accept(\Delta C, T): accept routine
# perturb(F): perturbation routine
def sa(Tmin, Tmax, N, alpha, F, costEval, accept, perturb, plot = True):
  assert(alpha < 1. and Tmin < Tmax)
  T = Tmax
  C = costEval(F)
  minC = C
  minF = F[:]
  Clist = []
  Temp = []
  while T > Tmin:
    for i in range(N):
      Fnew = perturb(F)
      Cnew = costEval(Fnew)
      if accept(Cnew - C, T):
        C, F = Cnew, Fnew
        if minC >= Cnew:
          minC, minF = Cnew, Fnew
        Clist.append(Cnew)
        Temp.append(T)
    T = T * alpha
  if plot:
    import matplotlib.pyplot as plt
    plt.plot(Temp, Clist)
    plt.xlim(max(Temp), min(Temp))
    plt.xscale('log')
  return minF, minC

def plot(coords):
  import matplotlib.pyplot as plt
  from matplotlib.patches import Rectangle

  fig, ax = plt.subplots()
  ax.plot([0, 0])
  ax.set_aspect('equal')
  ax.set_xlim(0,max([r[0][0] + r[1][0] for r in coords]))
  ax.set_ylim(0,max([r[0][1] + r[1][1] for r in coords]))

  for i,r in enumerate(coords):
    match i%4:
      case 3: hatch, color = '/+', 'red'
      case 2: hatch, color = '///', 'green'
      case 1: hatch, color = '/\//\//\/', 'blue'
      case _: hatch, color = '\\\\\\', 'gray'
    ax.add_patch(Rectangle(r[0], r[1][0], r[1][1],
              edgecolor = color, facecolor=color,
              hatch=hatch, fill = False,
              lw=2))
    ax.text(r[0][0] + r[1][0]//2, r[0][1] + r[1][1]//2, r[2], fontsize=8)
  plt.show()

class Module:
    def __init__(self, name, area, aspect_ratios):
      self._name = name
      self._area = area
      self._xy = [(math.sqrt(area*r), math.sqrt(area/r)) for r in aspect_ratios]
    def __repr__(self):
      return f"'{self._name} area:{self._area} xy:{self._xy}'"


### FLOORPLANNING USING SEQUENCE PAIRS ###
    
# m is a list of modules; see definition of the module class above
def sp_floorplan(m):

  def perturb_sp(Fi):
    F = copy.deepcopy(Fi)
    done = False
    while(done == False):
      choice = random.randint(1,3)
      last = len(F[0]) - 1
      match choice:
        case 1: # Swap blocks in F[0] : gamma+
          rand_ind = (0,0)
          while (rand_ind[0] == rand_ind[1]):
            rand_ind = (random.randint(0, last), random.randint(0, last))
          F[0][rand_ind[0]], F[0][rand_ind[1]] = F[0][rand_ind[1]], F[0][rand_ind[0]]
          done = True

        case 2: # Swap blocks in F[0] and F[1]: gamma+ and gamma-
          rand_ind = (0,0)
          while (rand_ind[0] == rand_ind[1]):
            rand_ind = (random.randint(0, last), random.randint(0, last))

          rand_ind2 = [0,0]
          for i in range(len(F[1])):
            if F[1][i]._name ==F[0][rand_ind[0]]._name: rand_ind2[0] = i
            if F[1][i]._name ==F[0][rand_ind[1]]._name: rand_ind2[1] = i
          F[0][rand_ind[0]], F[0][rand_ind[1]] = F[0][rand_ind[1]], F[0][rand_ind[0]]
          F[1][rand_ind2[0]], F[1][rand_ind2[1]] = F[1][rand_ind2[1]], F[1][rand_ind2[0]]
          done = True
          
        case _: # Change aspect ratio of a block
          rand_ind = random.randint(0, last)
          if len(F[0][rand_ind]._xy) > 1:
            rand_ind_xy = random.randint(0, len(F[0][rand_ind]._xy) - 1)
            F[0][rand_ind]._asp_r = rand_ind_xy
            done = True
    return F

  def costEval_sp(F):
    modules = findCoord_sp(F)
    x_max = 0
    y_max = 0
    for keys in modules:
      x_max = max(x_max, modules[keys]._x + modules[keys]._xy[modules[keys]._asp_r][0])
      y_max = max(y_max, modules[keys]._y + modules[keys]._xy[modules[keys]._asp_r][1])
    return x_max * y_max

  # class for storing the constraint graph table and block position
  class Mod_Pos:
    def __init__(self, name, area, xy):
      self._name = name
      self._area = area
      self._xy = xy
      self._left = []
      self._below = []
      self._x = 0
      self._y = 0
      self._asp_r = 0

  def findCoord_sp(F):
    modules = dict()
    for module in F[0][:]:
      modules[module._name] = module

    # horizontal & vertical constraint graph
    for module in F[0][:]:
      module._x = 0
      module._y = 0

    for module in F[0][:]:
      i,j = 0,0
      for k in range(len(F[1])):
        if F[0][k]._name ==module._name: i = k
        if F[1][k]._name ==module._name: j = k

      after_pos = set(F[0][i+1:])
      after_neg = set(F[1][j+1:])

      left = list(set(F[0][:]) - (after_pos | after_neg) - set([module]))
      below = list(after_pos & (set(F[1][:]) - after_neg))

      modules[module._name]._left = [modules[item._name] for item in left]
      modules[module._name]._below = [modules[item._name] for item in below]

    list_mod_h = F[0][:]
    list_mod_v = F[0][:]

    q_h = Queue(maxsize = len(F[0]))
    q_v = Queue(maxsize = len(F[1]))
    
    for iter in range(len(list_mod_h)):
      list_remove_h = []
      list_remove_v = []

      for i in list_mod_h:
        if modules[i._name]._left == []:
          q_h.put_nowait(modules[i._name])
          list_remove_h.append(i)

      for i in list_mod_v:
        if modules[i._name]._below == []:
          q_v.put_nowait(modules[i._name])
          list_remove_v.append(i)

      for i in list_remove_h:
        list_mod_h.remove(i)

      for i in list_remove_v:
        list_mod_v.remove(i)

      # Calculating positions of children blocks
      while q_h.qsize(): 
        deque_ele_h = q_h.get_nowait()
        for k in list_mod_h:
          if deque_ele_h in k._left:
            modules[k._name]._x = max(modules[k._name]._x, deque_ele_h._x + deque_ele_h._xy[deque_ele_h._asp_r][0])

      while q_v.qsize(): 
        deque_ele_v = q_v.get_nowait()
        for k in list_mod_v:
          if deque_ele_v in k._below:
            modules[k._name]._y = max(modules[k._name]._y, deque_ele_v._y + deque_ele_v._xy[deque_ele_v._asp_r][1])

      # Removing blocks from mod._left and mod._below once a block is dequed
      for i in list_remove_h:
        for mod in list_mod_h:
          if i in mod._left:
            mod._left.remove(i)

      for i in list_remove_v:
        for mod in list_mod_v:
          if i in mod._below:
            mod._below.remove(i)
    return modules

# MAIN
  modules = dict()
  for module in m:
    modules[module._name] = Mod_Pos(module._name, module._area, module._xy)
  ind = [modules[key] for key in modules]
  random.shuffle(ind)

  F = [ind[:],ind[:]]
  Tmax = sum([i._area for i in m])
  Fmin, Cmin = sa(1, Tmax, 10, 0.9, F, costEval_sp, accept, perturb_sp, False)
  modules = findCoord_sp(Fmin)
  coords = [((modules[k]._x,modules[k]._y),modules[k]._xy[modules[k]._asp_r],modules[k]._name) for k in modules]
  
  return (coords, Cmin)