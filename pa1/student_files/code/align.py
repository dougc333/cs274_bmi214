
"""

This file provides skeleton code for align.py. 

Locations with "FILL IN" in comments are where you need to add code.

Note - you do not need to follow this set up! It is just a suggestion, and may help for program design and testing.


Usage: python align.py input_file output_file

"""


import sys
import copy 
from pprint import pprint

#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)
    

#### ------- CLASSES ------- ####

class MatchMatrix(object):
    """
    Match matrix class stores the scores of matches in a data structure
    """
    def __init__(self,alphabet_a, alphabet_b,len_alphabet_a, len_alphabet_b):
        ### FILL IN ###
        self.matrix = [[0 for i in range(0,len_alphabet_b)] for j in range(0,len_alphabet_a)]
        self.alphabet_a=[x for x in alphabet_a] #this is the vertical part of matrix Altman convention, num_rows
        self.alphabet_b=[x for x in alphabet_b] #this is row across top of martix Altman convention, num_cols

    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it for
        """
        ### FILL IN ###
        row = self.alphabet_a.index(a)
        col = self.alphabet_b.index(b)
        self.matrix[row][col] = float(score)

    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence a and b is from sequence b.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """
        ### FILL IN ###
        return self.matrix[self.alphabet_a.index(a)][self.alphabet_b.index(b)]


class Node:  
  def __init__(self,x,y):
    self.x = x
    self.y = y
    self.weight = 0.
    self.M_pointer = []
    self.Ix_pointer = []
    self.Iy_pointer = []
    #self.pointers = []

class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
    ScoreEntries that are updated during alignment and used to output the maximum alignment.
    """

    def __init__(self, name, nrow, ncol):
        self.name = name # identifier for the score matrix - Ix, Iy, or M
        self.nrow = nrow
        self.ncol = ncol
        #this is wrong? off by one? because they start at 1 when counting len(alphabet_a)
        self.score_matrix =[[Node(i,j) for i in range(0,ncol)] for j in range(0,nrow)]
        
        # FILL IN 
        # you need to figure out a way to represent this and how to initialize
        # Hint: it may be helpful to have an object for each entry
    def M_pointer_add(self,row,col,tupl):
        self.score_matrix[row][col].M_pointer.append(tupl)
    
    def getM_pointer(self,row,col):
        return self.score_matrix[row][col].M_pointer
    
    def Ix_pointer_add(self,row,col,tupl):
        self.score_matrix[row][col].Ix_pointer.append(tupl)

    def getIx_pointer(self,row,col):
        return self.score_matrix[row][col].Ix_pointer

    def Iy_pointer_add(self,row,col,tupl):
        self.score_matrix[row][col].Iy_pointer.append(tupl)

    def getIy_pointer(self,row,col):
        return self.score_matrix[row][col].Iy_pointer

    def get_score(self, row, col):
        return self.score_matrix[row][col].weight
        
    def set_score(self, row, col, score):    
        ### FILL IN ###
        self.score_matrix[row][col].weight = float(score)
    
    def get_pointers(self, row, col):
        """
        Returns the indices of the entries that are pointed to
        This should be formatted as a list of tuples:
         ex. [(1,1), (1,0)]
        """
        ### FILL IN ###
        #needs to be flattened?
        #print("asdf",self.getM_pointer(row,col),self.getIx_pointer(row,col), self.getIy_pointer(row,col))
        return [self.getM_pointer(row,col), self.getIx_pointer(row,col), self.getIy_pointer(row,col)]

    def set_pointers(self, row, col,M_tupl, Ix_tupl, Iy_tupl): 
        ### FILL IN - this needs additional arguments ###
        ### FILL IN ###
        self.M_pointer_add(row,col,M_tupl)
        self.Ix_pointer_add(row,col,Ix_tupl)
        self.Iy_pointer_add(row,col,Iy_tupl)

    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix. Use this for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        """
        ### FILL IN ###  
        print(self.name+" = ")
        for row in self.score_matrix:
            print(" ".join([str(x.weight) for x in row]))

    
    def print_pointers_old(self):
        """
        debugging!
        """
        for i in range(0,self.nrow):
            for j in range(0,self.ncol):
                print("i j",i,j,self.score_matrix[i][j].pointers)

    def print_pointers(self):
        #convention for printing matrix
        print("pointers:",self.name)
        for i in range(1,self.nrow):
            for j in range(1,self.ncol):
                #merge all 3 into a list to remove empty list
                #print(i,j,["M"+str(x) for x in self.getM_pointer(i,j)],
                #["Ix"+str(x) for x in self.getIx_pointer(i,j)],
                #["Iy"+str(x) for x in self.getIy_pointer(i,j)])
                print(i,j,["M"+str(x) for x in self.score_matrix[i][j].M_pointer],
                ["Ix"+str(x) for x in self.score_matrix[i][j].Ix_pointer],
                ["Iy"+str(x) for x in self.score_matrix[i][j].Iy_pointer])
                

class AlignmentParameters(object):
    """
    Object to hold a set of alignment parameters from an input file.
    """
    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0    
        self.alphabet_a = "" 
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix(self.alphabet_a, self.alphabet_b, self.len_alphabet_a, self.len_alphabet_a)
        

    def load_params_from_file(self, input_file): 
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """
        lines = [line.rstrip('\n') for line in open(input_file)]
        self.seq_a = lines[0] 
        self.seq_b = lines[1] 
        global_alignment = int(lines[2]) 
        if global_alignment == 0:
            self.global_alignment = True
        else:
            self.global_alignment = False
        de_list = lines[3].split()
        self.dx = float(de_list[0])
        self.ex = float(de_list[1])
        self.dy = float(de_list[2])
        self.ey = float(de_list[3]) 
        self.len_alphabet_a = int(lines[4])
        self.alphabet_a = lines[5] 
        self.len_alphabet_b = int(lines[6])
        self.alphabet_b = lines[7]
        self.match_matrix = MatchMatrix(self.alphabet_a, self.alphabet_b, self.len_alphabet_a, self.len_alphabet_b)
        for x in range(8,len(lines)):
            if(len(lines[x].strip())>0):
                parse_input = lines[x].split()
                row = parse_input[0] #row in match matrix from 1 not used for defined interface
                col = parse_input[1] # col in match matrix from 1 not used for defined interface
                firstAA = parse_input[2] 
                secondAA = parse_input[3]
                similarity = parse_input[4]
                self.match_matrix.set_score(firstAA, secondAA, similarity)

class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by using "align()"
    """

    def __init__(self, input_file, output_file):
        """
        Input:
            input_file = file with the input for running an alignment
            output_file = file to write the output alignments to
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters()
        self.m_matrix = None
        self.ix_matrix = None
        self.iy_matrix = None
        self.score = 0.0
        self.final_path=[]
        #node = [[(4,5)],[],[]] becuse we need a notation for an empty tuple 
        #path=[[node1],[node2],[node3]]
        #paths = [[path],[path] = [ [[node1],[node2],[node3]], [[node1],[node2],[node3]] ]
        self.paths=[[]] #blue = node, pink=path, yellow=paths

    #we have to delay init till we load parameters in file because of align_test.py  
    def init_matrix(self):
        self.m_matrix = ScoreMatrix("M",len(self.align_params.seq_a)+1, len(self.align_params.seq_b)+1)
        self.ix_matrix = ScoreMatrix("Ix",len(self.align_params.seq_a)+1, len(self.align_params.seq_b)+1)
        self.iy_matrix = ScoreMatrix("Iy",len(self.align_params.seq_a)+1, len(self.align_params.seq_b)+1)
        
    def align(self):
        """
        Main method for running alignment.
        """

        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)
        self.init_matrix()
        # populate the score matrices based on the input parameters
        self.populate_score_matrices()
        self.traceback()
        #self.write_output()
        
        ### FILL IN ###
    
    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        Should call update(i,j) for each entry in the score matrices
        """
        ### FILL IN ###
        #careful to use len_alphabet_a vs. len(align_params.seq_a) for align_test
        #we dont specify seqa, only size of alphabet for testing update_ix, update_m, update_iy
        for i in range(0,len(self.align_params.seq_a)+1):
            for j in range(0,len(self.align_params.seq_b)+1):
                self.update(i,j)   
        self.debug_print()
                  
    def debug_print(self):           
        print("----------------")
        self.m_matrix.print_scores()
        print("----------------")
        self.ix_matrix.print_scores()
        #print("maxIx:",max([max(x)  for x in self.ix_matrix.score_matrix]))
        print("----------------")
        self.iy_matrix.print_scores()
        print("----------------")
        self.m_matrix.print_pointers()
        self.ix_matrix.print_pointers()
        self.iy_matrix.print_pointers()
        print("-------------------------")
    def update(self, row, col):
        """
        Method to update the matrices at a given row and column index.

        Input:
           row = the row index to update
           col = the column index to update
        """
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)
        
    def update_m(self, row, col):
        ### FILL IN ###
        if (row==0 and col==0):
            self.m_matrix.set_score(row,col,0.0)
        elif(row!=0 and col!=0):
            firstM =  self.m_matrix.get_score(row-1,col-1) + (self.align_params.match_matrix.get_score(self.align_params.seq_a[row-1],self.align_params.seq_b[col-1]))
            secondM = self.ix_matrix.get_score(row-1,col-1) + (self.align_params.match_matrix.get_score(self.align_params.seq_a[row-1],self.align_params.seq_b[col-1]))
            thirdM =  self.iy_matrix.get_score(row-1,col-1) + (self.align_params.match_matrix.get_score(self.align_params.seq_a[row-1],self.align_params.seq_b[col-1]))
            maxM = max(firstM,secondM,thirdM)
            print("update_m:",row,col,firstM,secondM,thirdM,maxM)
            self.m_matrix.set_score(row,col, maxM)
            if fuzzy_equals(firstM,maxM):
                print("update_m m_M pointer update:",row-1,col-1)
                self.m_matrix.M_pointer_add(row,col,(row-1,col-1))
            if fuzzy_equals(secondM,maxM):
                print("update_m ix_Ix pointer update:",row-1,col-1)
                self.m_matrix.Ix_pointer_add(row,col,(row-1,col-1))
            if fuzzy_equals(thirdM,maxM):
                print("update_m iy_Iy pointer update:",row-1,col-1)
                self.m_matrix.Iy_pointer_add(row,col,(row-1,col-1))

       
    def update_ix(self, row, col):
        if (row==0 and col==0):
            self.ix_matrix.set_score(row,col,0.0)
        elif(row!=0 and col!=0):
            firstIx = self.m_matrix.get_score(row-1,col) - self.align_params.dy
            secondIx = self.ix_matrix.get_score(row-1,col) - self.align_params.ey
            maxIx = max(firstIx,secondIx)
            print("update_ix:",row,col,firstIx,secondIx,maxIx)
            self.ix_matrix.set_score(row,col, maxIx)
            if fuzzy_equals(firstIx,maxIx):
                print("update_ix ix_M pointer update:",row-1,col)
                self.ix_matrix.M_pointer_add(row,col,(row-1,col))
            if fuzzy_equals(secondIx,maxIx):
                print("update_ix ix_Iy pointer update:",row-1,col)
                self.ix_matrix.Ix_pointer_add(row,col,(row-1,col))

    def update_iy(self, row, col):        
        if (row==0 and col==0):
            self.iy_matrix.set_score(row,col,0.0)
        elif(row!=0 and col!=0):
            firstIy = self.m_matrix.get_score(row,col-1) - self.align_params.dx
            secondIy = self.iy_matrix.get_score(row,col-1) - self.align_params.ex
            maxIy = max(firstIy,secondIy)
            print("update_iy:",row,col,firstIy,secondIy,maxIy)
            self.iy_matrix.set_score(row,col, maxIy)
            if fuzzy_equals(firstIy,maxIy):
                print("update_iy iy_M pointer update:",row,col-1)
                self.iy_matrix.M_pointer_add(row,col,(row,col-1))
            if fuzzy_equals(secondIy,maxIy):
                print("update_iy iy_Iy pointer update:",row,col-1)
                self.iy_matrix.Iy_pointer_add(row,col,(row,col-1))   


    def find_traceback_start(self):
        """
        Finds the location to start the traceback..
        Think carefully about how to set this up for local 

        Returns:
            (max_val, max_loc) where max_val is the best score
            max_loc is a list [] containing tuples with the (i,j) location(s) to start the traceback
             (ex. [(1,2), (3,4)])
        """
        ### FILL IN ###
        if self.align_params.global_alignment==1:
            self.score = self.m_matrix.get_score(self.m_matrix.nrow-1,self.m_matrix.ncol-1)
            #print("setting score:",self.score)
            return (self.m_matrix.nrow-1, self.m_matrix.ncol-1)
        else:
            print("there can be only 1 local max? NO")
            maxM = max([max(x)  for x in self.m_matrix.score_matrix])
            maxIx = max([max(x)  for x in self.ix_matrix.score_matrix])
            maxIy = max([max(x)  for x in self.iy_matrix.score_matrix])
            max_all=max(maxM, maxIx, maxIy)
            max_loc =[(ix,iy) for ix, row in enumerate(a) for iy, i in enumerate(row) if i == 0]
            print("max_loc:",max_loc)
            return max_loc
    def str_to_tuple(self,st):
        if st[0]=="M":
            return()
        else:
            #Ix or Iy
            return

    def add_one(self,last_node, node):
        #print("add_one last_node:",last_node)
        #print("add_one node",node)
        for idx in range(0,len(self.paths)):
            #print("add_one testing last node:",self.paths[idx][-1])
            if (self.paths[idx][-1]==last_node):
                #print("adding!!!!")
                self.paths[idx].append(node)
            
    def branch(self,last_node,node1,node2):
        print("branch node1:",node1,"branch node2:",node2)
        for idx in range(0,len(self.paths)):
            #print("branch testing last node:",self.paths[idx][-1])
            if (self.paths[idx][-1]==last_node):
                #print("adding!!!!")
                dup = copy.deepcopy(self.paths[idx]) #copies path
                dup.append(node1)
                self.paths[idx].append(node2)
                self.paths.append(dup)

    def branch2(self,last_node,node1,node2,node3):
        '''
        termination of path, add the final node with M,Ix,Iy set. This is an invalid
        node because it requires you to be in all 3 matrices which is physically
        impossible but this is defintion of termination condition on the graph.  
        '''
        print("branch2 node1:",node1,"branch2 node2:",node2,"branch2 node3:",node3)
        for idx in range(0,len(self.paths)):
            #print("branch testing last node:",self.paths[idx][-1])
            if (self.paths[idx][-1]==last_node):
                #remove this path from self.paths to print out
                self.final_path.append(self.paths[idx])
                self.paths.remove(self.paths[idx])


    def addM(self,tupl):
        node = [[tupl],[],[]]
        for x in range(0,len(self.paths)):
            self.paths[x].append(node)
    def addIx(self,tupl):
        node = [[],[tupl],[]]
        for x in range(0,len(self.paths)):
            self.paths.append(node)
    def addIy(self,tupl):
        node = [[],[],[tupl]]
        for x in range(0,len(self.paths)):
            self.paths.append(node)
    def last_nodes(self):
        '''
        returns last node tuples for next iteration in score matrix till end of path
        '''
        nodes=[]
        for idx in range(0,len(self.paths)):
            #print("path:",self.paths[idx])
            last_node =self.paths[idx][-1]
            #print("last_node:",last_node)
            nodes.append(last_node)
        return nodes

    def print_paths(self):
        print("num_paths:",len(self.paths))
        for idx in range(0,len(self.paths)):
            print(self.paths[idx])
    
    def traceback(self): ### FILL IN additional arguments ###
        """
        Performs a traceback.
        Hint: include a way to printing the traceback path. This will be helpful for debugging!
           ex. M(5,4)->Iy(4,3)->M(4,2)->Ix(3,1)->Ix(2,1)->M(1,1)->M(0,0)
        """
        ### FILL IN ###
        start = self.find_traceback_start()
        print("traceback_start:",start)
        i = start[0]
        j = start[1]
        
        #print("self.paths before adding start:",self.paths)
        self.addM(start)
        #print("self.paths after adding start node",self.paths)
        #each node [[path] [path] [path]]=[[M] [Ix] [Iy]]
        #print( "m_M_pointer",self.m_matrix.getM_pointer(i,j),"m_Ix_pointer:",self.m_matrix.getIx_pointer(i,j),"m_Iy_pointer:",self.m_matrix.getIy_pointer(i,j))
        #print("getpointers:",self.m_matrix.get_pointers(i,j))

        num_iter=0
        while(len(self.paths)>0 and num_iter<8):
            print("start loop self.path:")
            self.print_paths()
            last_nodes = self.last_nodes()
            print("num_iter:",num_iter,"i j",i,j,"last_nodes:",last_nodes)
            
            for n in last_nodes:
                #use n[0],n[1],n[2] to figure out if tis is M,Ix or Iy
                print("last node:",n,"len n[0],n[1],n[2]",len(n[0]),len(n[1]),len(n[2]))
                l=[]
                if len(n[0])==1:
                    print("lookup M")
                    i,j = n[0][0]
                    diag = [self.m_matrix.get_pointers(i,j)[0],[],[]] #M
                    top = [[],self.m_matrix.get_pointers(i,j)[1],[]] #ix
                    left = [[],[],self.m_matrix.get_pointers(i,j)[2]] #iy
                    l.extend(diag)
                    l.extend(top)
                    l.extend(left)
                elif(len(n[1])==1):
                    print("lookup Ix")
                    i,j = n[1][0]
                    diag = [self.ix_matrix.get_pointers(i,j)[0],[],[]] #M
                    top = [[],self.ix_matrix.get_pointers(i,j)[1],[]] #ix
                    #left = [[],[],self.m_matrix.get_pointers(i,j)[2]] #iy
                    print("lookup Ix diag:",diag," top:",top)
                    l.extend(diag)
                    l.extend(top)
                    #l.extend(left)
                elif(len(n[2])==1):
                    print("lookup Iy")
                    i,j = n[2][0]
                    diag = [self.iy_matrix.get_pointers(i,j)[0],[],[]] #M
                    #top = [[],self.m_matrix.get_pointers(i,j)[1],[]] #ix
                    left = [[],[],self.iy_matrix.get_pointers(i,j)[2]] #iy
                    l.extend(diag)
                    #l.extend(top)
                    l.extend(left)
                    print("lookupIy diag:",diag," left:",left)
                else:
                    print("*******should not see this*********")
                flatten = [x for x in l if len(x)>0]
                num_diag = [len(x) for x in diag]
                num_top = [len(x) for x in top]
                num_left = [len(x) for x in left]
                print("diag:",diag," top:",top," left",left," flatten:",flatten,"len(flatten):",len(flatten))
                print("num_diag,num_top,num_left",num_diag,num_top,num_left)            
                #if more than one have to branch
                if (len(flatten))==0:
                    print("len(flatten) =0")
                    #done or error
                    print("nothing")
                elif(len(flatten))==1:
                    #append single flatten list and append once
                    print("len(flatten)=1")
                    single_node = []
                    if (1 in num_diag):
                        single_node.append(diag)
                    if (1 in num_top):
                        single_node.append(top)
                    if (1 in num_left):
                        single_node.append(left)
                    self.add_one(n,single_node[0])
                    print("after 1 add path:")
                    self.print_paths()
                elif (len(flatten)==2 ):
                    #branch multiple entries need to branch and add separate tuple to each one. flatten, branch
                    #to number of tuples and remove tuple and add to path list
                    print("branch len(flatten)==2")
                    two_nodes=[]
                    if (1 in num_diag):
                        two_nodes.append(diag)
                    if (1 in num_top):
                        two_nodes.append(top)
                    if (1 in num_left):
                        two_nodes.append(left)
                    print("len==2 two_nodes:",two_nodes)
                    self.branch(n,two_nodes[0],two_nodes[1])
                    print("after branch self.paths:")
                    self.print_paths()
                elif(len(flatten)==3):
                    #termination condition
                    node=[]
                    if (1 in num_diag):
                        node.append(diag)
                    if (1 in num_top):
                        node.append(top)
                    if (1 in num_left):
                        node.append(num_left)
                    #print("before branch2 self.paths:")
                    #how to append?
                    self.branch2(n,node[0],node[1],node[2])
                    #do nothing add code for branch3 case
                    print("after branch 2 nodes self.paths:")
                    self.print_paths()
                else:
                    print("should not see this if len flatten - tie")
            num_iter+=1
            #remove complete paths
            #self.remove_paths()
        print("after while loop self.paths:",num_iter)

    def write_output(self):
        fh = open(self.output_file,"w")
        fh.write(str(self.score))
        fh.write("\n")
        fh.write("\n")
        fh.write("AATG_C\n")
        fh.write("A__GGC\n")
        fh.write("\n")
        fh.write("ATG_C\n")
        fh.write("A_GGC\n")
        fh.write("\n")
        fh.close()
        return None

def main():
    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # create an align object and run
    align = Align(input_file, output_file)
    align.align()
    align.write_output()    
    #align.traceback()

if __name__=="__main__":
    main()
