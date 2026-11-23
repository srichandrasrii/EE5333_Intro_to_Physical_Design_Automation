import numpy as np
class kdtree:
  def __init__(self, rects):
    self._rects = rects

  def region_query(self, region):
    ll_reg = region.ll
    ur_reg = region.ur
    nbrs = self._rects # actual nbrs

    while(len(nbrs) > 20):

      ll_x_median = np.median([rect.ll.x for rect in nbrs])
      ll_y_median = np.median([rect.ll.y for rect in nbrs])
      ur_x_median = np.median([rect.ur.x for rect in nbrs])
      ur_y_median = np.median([rect.ur.y for rect in nbrs])

      if(ur_reg.x < ll_x_median):
        sub_rects = [rect for rect in nbrs if rect.ll.x < ll_x_median]
      elif(ur_reg.y < ll_y_median):
        sub_rects = [rect for rect in nbrs if rect.ll.y < ll_y_median]
      elif(ll_reg.x > ur_x_median):
        sub_rects = [rect for rect in nbrs if rect.ur.x > ur_x_median]
      elif(ll_reg.y > ur_y_median):
        sub_rects = [rect for rect in nbrs if rect.ur.y > ur_y_median]
      else:
        break

      nbrs = sub_rects

    return nbrs
