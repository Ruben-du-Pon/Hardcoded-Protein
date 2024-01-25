from ..classes.protein import Protein
from ..classes.protein import Aminoacid
from .bfs import BfsFold
import random
from typing import List
import time
import numpy as np


class MctsFold(BfsFold):
    def __init__(self, protein: Protein, dimensions: int, when_cutting=6, step=1):
        """
        Initialize Mcts instance.

        Parameters:
        - protein (Protein): The protein structure to be folded.
        - dimensions (int): Folding done in 2D or 3D.
        - when_cutting (int): The length at which to start cutting the protein sequence during folding.
        - step (int): The step size to use during folding.
        """
        super().__init__(protein, dimensions, when_cutting, step)
        self._min_keys = []



    def __bfsfold(self, protein: Protein, when_cutting, step) -> List[str]:
        """
        Perform Breadth-First Search (BFS) based folding on the given protein structure.

        Parameters:
        - protein (Protein): The protein structure to be folded.
        - when_cutting (int): The length at which to start cutting the protein sequence during folding.
        - step (int): The step size to use during folding.
        - dimension (int): The dimension of the protein folding (default is 2).

        Returns:
        - Protein: The folded protein structure.
        """
        posit = [(0, 0, 0)]
        sequence_protein = protein._sequence
        if self.dimensions == 2:
            types = {"R", "L", "U", "D"}
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0)}
        elif self.dimensions == 3:
            types = {"R", "L", "U", "D", "F", "B"}
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0), "F": (0, 0, 1), "B": (0, 0, -1)}
        min_keys = set()
        
        if "C" in self._sequence and self.dimensions == 2:
            when_cutting = 6
        elif self.dimensions == 3:
            when_cutting = 5
        else:
            when_cutting = 5

        while len(self._sequence) <= when_cutting:
                when_cutting -= 1
        
        step = 1

        if len(self._protein) <= 10:
            going_till = len(self._protein) - 1
        else:
            going_till = 10

        for depth in range(when_cutting, going_till, step):
            create_d = super()._BfsFold__create_dict(
                protein, sequence_protein, types, depth, step, min_keys, posit
            )
            posit = [(0, 0, 0)]

            min_key = min(create_d, key=lambda k: create_d[k])
            min_keys = {k for k, v in create_d.items() if v == create_d[min_key]}
            if len(min_keys) >= 2:
                # Randomly select 1 of the set
                min_keys = random.sample(min_keys, 1)
            else:
                min_keys = list(min_keys)
            for key in min_keys[0]:
                posit.append(tuple(np.array(posit[-1]) + np.array(move[key])))


        return min_keys
    

    def __get_protein_score(self, min_keys):
        min_keys_ = min_keys[0]

        if self.dimensions == 2:
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0)}
        elif self.dimensions == 3:
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0), "F": (0, 0, 1), "B": (0, 0, -1)}

        prt = Protein(self._sequence)
        current = prt.get_head()
        current.position = (0, 0, 0)

        prt.add_to_grid(current.position, current)

        for key in min_keys_:
            current = current.link
            current.position = tuple(np.array(current.predecessor.position) + np.array(move[key]))
            prt.add_to_grid(current.position, current)

        return prt.get_score()


    def __create_final_protein(self, min_keys):
        min_keys_ = min_keys[0]

        if self.dimensions == 2:
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0)}
        elif self.dimensions == 3:
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0), "F": (0, 0, 1), "B": (0, 0, -1)}

        prt = Protein(self._sequence)
        current = prt.get_head()

        prt.add_to_grid(current.position, current)

        for key in min_keys_:
            current = current.link
            current.position = tuple(np.array(current.predecessor.position) + np.array(move[key]))
            prt.add_to_grid(current.position, current)

        return prt

    
    def __get_coordinates(self, min_keys):
        # print(min_keys)
        pos = [(0, 0, 0)]
        if self.dimensions == 2:
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0)}
        elif self.dimensions == 3:
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0), "F": (0, 0, 1), "B": (0, 0, -1)}

        for key in min_keys[0]:
            pos.append(tuple(np.array(pos[-1]) + np.array(move[key])))

        return pos
    


    def __mcts(self, min_keys):
        length_protein, min_keys_, dict_scores = len(self._protein), min_keys, {}

        if self.dimensions == 2:
            types = {"R", "L", "U", "D"}
            types_ = {"R", "L", "U", "D"}
        elif self.dimensions == 3:
            types = {"R", "L", "U", "D", "F", "B"}
            types_ = {"R", "L", "U", "D", "F", "B"}

        while length_protein != (len(min_keys[0]) + 1):
            if min_keys_[0][-1] == "R":
                types.remove("L")
            elif min_keys_[0][-1] == "L":
                types.remove("R")
            elif min_keys_[0][-1] == "U":
                types.remove("D")
            elif min_keys_[0][-1] == "D":
                types.remove("U")
            elif min_keys_[0][-1] == "F":
                types.remove("B")
            elif min_keys_[0][-1] == "B":
                types.remove("F")

            for type in types:
                for iteration in range(100000):
                    while length_protein != (len(min_keys_[0]) + 1):
                        if self.dimensions == 2:
                            types_ = {"R", "L", "U", "D"}
                        elif self.dimensions == 3:
                            types_ = {"R", "L", "U", "D", "F", "B"}

                        if min_keys_[0][-1] == "R":
                            types_.remove("L")
                        elif min_keys_[0][-1] == "L":
                            types_.remove("R")
                        elif min_keys_[0][-1] == "U":
                            types_.remove("D")
                        elif min_keys_[0][-1] == "D":
                            types_.remove("U")
                        elif min_keys_[0][-1] == "F":
                            types_.remove("B")
                        elif min_keys_[0][-1] == "B":
                            types_.remove("F")

                        min_keys_ = [min_keys_[0] + random.choice(list(types_))]
                    
                    if iteration == 0:
                        dict_scores[type] = self.__get_protein_score(min_keys_)
                    else:
                        dict_scores[type] = dict_scores[type] + self.__get_protein_score(min_keys_)

                    min_keys_ = min_keys
                
                dict_scores[type] = dict_scores[type] / 50000
                print(dict_scores)
            
            min_keys__ = [min_keys[0] + min(dict_scores, key=dict_scores.get)]
            coordinates_ = self.__get_coordinates(min_keys__)

            while len(coordinates_) != len(set(coordinates_)):
                # print(1, coordinates_)
                del dict_scores[min(dict_scores, key=dict_scores.get)]
                # print(2, dict_scores)
                min_keys__ = [min_keys[0] + min(dict_scores, key=dict_scores.get)]
                coordinates_ = self.__get_coordinates(min_keys__)

            min_keys = min_keys__
            dict_scores = {}

            if self.dimensions == 2:
                types = {"R", "L", "U", "D"}
                types_ = {"R", "L", "U", "D"}
            elif self.dimensions == 3:
                types = {"R", "L", "U", "D", "F", "B"}
                types_ = {"R", "L", "U", "D", "F", "B"}


        return self.__create_final_protein(min_keys)


    def run(self) -> Protein:
        start_time = time.time()

        if len(self._protein) > 10:
            min_keys = self.__bfsfold(self._protein, self._cut, self._step) 
            # print(min_keys)
            result = self.__mcts(min_keys)
        else:
            result = super()._BfsFold__bfsfold(self._protein, self._cut, self._step)
        end_time = time.time()  # Record the end time
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time} seconds")

        return result
                

                
