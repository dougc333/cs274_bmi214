
class AlignGraph:

  def __init__(self):
    self.first_string=[x for x in "TGTTA"]
    self.second_string = [x for x in "TCGT"]
    self.num_cols = len(self.first_string)
    self.num_rows = len(self.second_string)
    self.graph = []
    for i in range(0,num_rows)
      rows=[]
      for j in range(0,num_cols):
        rows.append(Node(i,j))
      self.graph.append(rows)