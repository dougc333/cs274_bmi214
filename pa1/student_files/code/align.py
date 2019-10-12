
"""

This file provides skeleton code for align.py. 

Locations with "FILL IN" in comments are where you need to add code.

Note - you do not need to follow this set up! It is just a suggestion, and may help for program design and testing.


Usage: python align.py input_file output_file

"""


import sys


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
        self.matrix[row][col] =score

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
    self.east=None
    self.se=None
    self.south=None
    self.weight = 0.
    self.pointers = []

  def set_weights(self,east, se, south):
    self.east = east
    self.se = se
    self.south = south

  def node_print():
    print("node_print:",self.x, self.y, self.east, self.se, self.south. self.weight)
    print("pointers:",pointers)


class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
    ScoreEntries that are updated during alignment and used to output the maximum alignment.
    """

    def __init__(self, name, nrow, ncol):
        self.name = name # identifier for the score matrix - Ix, Iy, or M
        self.nrow = nrow
        self.ncol = ncol
        self.score_matrix =[]
        for i in range(0,self.nrow):
            row=[]
            for j in range(0,self.ncol):
                row.append(Node(i,j))
            self.score_matrix.append(row)
        # FILL IN 
        # you need to figure out a way to represent this and how to initialize
        # Hint: it may be helpful to have an object for each entry

    def get_score(self, row, col):
        return self.score_matrix[row][col].weight
        
    def set_score(self, row, col, score):    
        ### FILL IN ###
        self.score_matrix[row][col].weight = score
    
    def get_pointers(self, row, col):
        """
        Returns the indices of the entries that are pointed to
        This should be formatted as a list of tuples:
         ex. [(1,1), (1,0)]
        """
        ### FILL IN ###
        return self.score_matrix[row][col].pointers

    def set_pointers(self, row, col,add_tuple): ### FILL IN - this needs additional arguments ###
        ### FILL IN ###
        self.score_matrix[row][col].pointers.append(tuple)

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
        #add node property printout
        #for i in range(0,self.nrow):
        #    for j in range(0,self.ncol):
        #        print("i j",i,j,self.score_matrix[i][j].weight)
        # 
        print(self.name+" = ")
        for row in self.score_matrix:
            print(" ".join([str(x.weight) for x in row]))

    
    def print_pointers(self):
        """
        Returns a nicely formatted string containing the pointers for each entry in the score matrix. Use this for debugging!
        """

        ### FILL IN ###
        for i in range(0,nrow):
            for j in range(0,ncol):
                print(self.score_matrix.)
    

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
        self.global_alignment = int(lines[2]) 
        de_list = lines[3].split()
        self.dx = float(de_list[0])
        self.ex = float(de_list[1])
        self.dy = float(de_list[2])
        self.ey = float(de_list[3]) 
        self.len_alphabet_a = int(lines[4])
        self.alphabet_a = lines[5] 
        self.len_alphabet_b = int(lines[6])
        self.alphabet_b = lines[7]
        #print("num_lines:",len(lines))
        self.match_matrix = MatchMatrix(self.alphabet_a, self.alphabet_b, self.len_alphabet_a, self.len_alphabet_b)
        for x in range(8,len(lines)):
            #print(lines[x])
            parse_input = lines[x].split()
            row = parse_input[0] #redundant
            col = parse_input[1] #redundant
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
        self.align_params.load_params_from_file(self.input_file) 

        ### FILL IN - note be careful about how you initialize these! ###
        self.M = ScoreMatrix("M",len(self.align_params.seq_a)+1, len(self.align_params.seq_b)+1)
        self.Ix = ScoreMatrix("Ix",len(self.align_params.seq_a)+1, len(self.align_params.seq_b)+1)
        self.Iy = ScoreMatrix("Iy",len(self.align_params.seq_a)+1, len(self.align_params.seq_b)+1)
        

    def align(self):
        """
        Main method for running alignment.
        """

        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()
        self.traceback()
        self.write_output()
        # perform a traceback and write the output to an output file

        ### FILL IN ###
    
    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        Should call update(i,j) for each entry in the score matrices
        """
        ### FILL IN ###
        for i in range(0,len(self.align_params.seq_a)+1):
            for j in range(0,len(self.align_params.seq_b)+1):
                if (i==0 and j==0):
                    self.M.score_matrix[i][j].weight=0.0
                    self.Ix.score_matrix[i][j].weight=0.0
                    self.Iy.score_matrix[i][j].weight=0.0
                elif(i==0 and j>0):
                    self.M.score_matrix[i][j].weight=0.0
                    self.Iy.score_matrix[i][j].weight = 0.0
                    self.Ix.score_matrix[i][j].weight = 0.0
                elif(j==0 and i>0):
                    self.M.score_matrix[i][j].weight=0.0
                    self.Iy.score_matrix[i][j].weight = 0.0
                    self.Ix.score_matrix[i][j].weight = 0.0
                elif(i!=0 and j!=0):
                    #updateM, updateIx, updateIy
                    #print(self.align_params.seq_a[i-1])
                    #print(self.align_params.seq_b[j-1])
                    #print(self.align_params.match_matrix.get_score(self.align_params.seq_a[i-1],self.align_params.seq_b[j-1]))
                    
                    self.M.score_matrix[i][j].weight = max(
                    self.M.score_matrix[i-1][j-1].weight + float(self.align_params.match_matrix.get_score(self.align_params.seq_a[i-1],self.align_params.seq_b[j-1])),
                    self.Ix.score_matrix[i-1][j-1].weight + float(self.align_params.match_matrix.get_score(self.align_params.seq_a[i-1],self.align_params.seq_b[j-1])),
                    self.Iy.score_matrix[i-1][j-1].weight + float(self.align_params.match_matrix.get_score(self.align_params.seq_a[i-1],self.align_params.seq_b[j-1]))
                    ) 
                    
                    self.Ix.score_matrix[i][j].weight  = max(
                    self.M.score_matrix[i-1][j].weight - self.align_params.dy,
                    self.Ix.score_matrix[i-1][j].weight - self.align_params.ey
                    )
                    self.Iy.score_matrix[i][j].weight = max(
                        self.M.score_matrix[i][j-1].weight - self.align_params.dx,
                        self.Iy.score_matrix[i][j-1].weight - self.align_params.ex
                    )
                    #set pointers to highest values between M, Ix,Iy as tuples of location ("M",x,y)
                    max_value = max(
                    self.M.score_matrix[i-1][j-1].weight,
                    self.Ix.score_matrix[i-1][j].weight,
                    self.Iy.score_matrix[i][j-1].weight 
                    )
                    
                    if max_value == self.M.score_matrix[i-1][j-1].weight:
                        self.M.score_matrix[i][j].pointers.append(("M(",i-1,j-1,")"))
                    elif(max_value == self.Ix.score_matrix[i-1][j].weight):
                        self.Ix.score_matrix[i][j].pointers.append(("Ix(",i-1,j,")"))
                    elif(max_value == self.Iy.score_matrix[i][j-1].weight):
                        self.Iy.score_matrix[i][j].pointers.append(("Iy(",i,j-1,")"))
                    else:
                        print("pointer setting should not see this i j",i,j)
                else:
                    print("should not see this")
        print("M")
        self.M.print_scores()
        print("----------------")
        print("Ix")
        self.Ix.print_scores()
        print("----------------")
        print("Iy")
        self.Iy.print_scores()
        print("----------------")
        self.M.score_matrix.print_pointers()
        #print("----------------")
        #self..score_matrix.print_pointers()
        #print("----------------")
        #self.M.score_matrix.print_pointers()
        
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
         self.M.score_matrix[i][j].weight = max(
            self.M.score_matrix[i-1][j-1].weight + float(self.align_params.match_matrix.get_score(self.align_params.seq_a[i-1],self.align_params.seq_b[j-1])),
            self.Ix.score_matrix[i-1][j-1].weight + float(self.align_params.match_matrix.get_score(self.align_params.seq_a[i-1],self.align_params.seq_b[j-1])),
            self.Iy.score_matrix[i-1][j-1].weight + float(self.align_params.match_matrix.get_score(self.align_params.seq_a[i-1],self.align_params.seq_b[j-1]))
        ) 
                    

    def update_ix(self, row, col):
        self.Ix.score_matrix[i][j].weight = max(
            self.M.score_matrix[i-1][j].weight - self.align_params.dy,
            self.Ix.score_matrix[i-1][j].weight - self.align_params.ey
            )

    def update_iy(self, row, col):
        self.Iy.score_matrix[i][j].weight = max(
            self.M.score_matrix[i][j-1].weight - self.align_params.dx,
            self.Iy.score_matrix[i][j-1].weight - self.align_params.ex
            )

    
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
            return self.align_params.len_alphabet_a+1, self.align_params.len_alphabet_b+1
        else:
            print("there can be only 1 local max?")
            

    def traceback(self): ### FILL IN additional arguments ###
        """
        Performs a traceback.
        Hint: include a way to printing the traceback path. This will be helpful for debugging!
           ex. M(5,4)->Iy(4,3)->M(4,2)->Ix(3,1)->Ix(2,1)->M(1,1)->M(0,0)
        """
        ### FILL IN ###
        return None

    def write_output(self):
        ### FILL IN ###
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


if __name__=="__main__":
    main()
