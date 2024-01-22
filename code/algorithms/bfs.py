from ..classes.protein import Protein
import random
from typing import List
import time


class BfsFold:
    def __init__(self, protein: Protein, dimensions: int, when_cutting=6, step=1):
        """
        Initialize BfsFold instance.

        Parameters:
        - protein (Protein): The protein structure to be folded.
        - dimensions (int): Folding done in 2D or 3D.
        - when_cutting (int): The length at which to start cutting the protein sequence during folding.
        - step (int): The step size to use during folding.
        """
        self._protein = protein
        self._sequence = protein._sequence
        self._cut = when_cutting
        self._step = step
        self.dimensions = dimensions


    def cutting_seq(self, protein: Protein, cut=10) -> List[str]:
        length, seq, proteins = len(protein), protein._sequence, []

        for x in range(0, length, 10):
            proteins.append(seq[x:x+cut]) 

        return proteins

    def __create_nested_dict(
        self, protein: Protein, keys, depth, prev=None, pos=[(0, 0, 0)]
    ):
        """
        Recursively create a nested dictionary representing possible folding positions.

        Parameters:
        - protein (Protein): The protein structure being folded.
        - keys (list): The possible folding directions (e.g., ["R", "L", "U", "D"]).
        - depth (int): The current depth in the folding process.
        - prev (str): The direction of the previous fold step.
        - pos (list): The list of positions representing the folded structure.

        Returns:
        - dict: A nested dictionary representing possible folding positions.
        """
        if self.dimensions == 2:
            aminoacid = protein._head

            if aminoacid.link is None or depth == 1:
                return {"pos": pos[-1]}

            result_dict = {"pos": pos[-1]}
            for key in keys:
                if (
                    (key == "R" and prev == "L")
                    or (key == "L" and prev == "R")
                    or (key == "U" and prev == "D")
                    or (key == "D" and prev == "U")
                ):
                    continue
                else:
                    move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0)}
                    result_dict[key] = self.__create_nested_dict(
                        protein,
                        keys,
                        depth - 1,
                        key,
                        pos + [tuple(x + y for x, y in zip(pos[-1], move[key]))],
                    )

            return result_dict
        
        elif self.dimensions == 3:
            aminoacid = protein._head

            if aminoacid.link is None or depth == 1:
                return {"pos": pos[-1]}

            result_dict = {"pos": pos[-1]}
            for key in keys:
                if (
                    (key == "R" and prev == "L")
                    or (key == "L" and prev == "R")
                    or (key == "U" and prev == "D")
                    or (key == "D" and prev == "U")
                    or (key == "F" and prev == "B")
                    or (key == "B" and prev == "F")
                ):
                    continue
                else:
                    move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0), "F": (0, 0, 1), "B": (0, 0, -1)}
                    result_dict[key] = self.__create_nested_dict(
                        protein,
                        keys,
                        depth - 1,
                        key,
                        pos + [tuple(x + y for x, y in zip(pos[-1], move[key]))],
                    )

            return result_dict 

    def __valid_combinations(self, keys, prev=None, length=2, it=0):
        """
        Generate valid combinations of folding directions.

        Parameters:
        - keys (list): The possible folding directions (e.g., ["R", "L", "U", "D"]).
        - prev (str): The direction of the previous fold step.
        - length (int): The desired length of the folding sequence.
        - it (int): Current iteration index.

        Returns:
        - list: List of valid combinations of folding directions.
        """
        if self.dimensions == 2:
            if it == length:
                return ['']

            valid_combos = []

            for key in keys:
                if (
                    (key == "R" and prev == "L")
                    or (key == "L" and prev == "R")
                    or (key == "U" and prev == "D")
                    or (key == "D" and prev == "U")
                ):
                    continue
                else:
                    combos = self.__valid_combinations(
                        keys, prev=key, length=length, it=it + 1
                    )
                    valid_combos.extend([key + combo for combo in combos])

            return valid_combos

        elif self.dimensions == 3:
            if it == length:
                return ['']

            valid_combos = []

            for key in keys:
                
                if (
                    (key == "R" and prev == "L")
                    or (key == "L" and prev == "R")
                    or (key == "U" and prev == "D")
                    or (key == "D" and prev == "U")
                    or (key == "F" and prev == "B")
                    or (key == "B" and prev == "F")
                ):
                    continue
                else:
                    combos = self.__valid_combinations(
                        keys, prev=key, length=length, it=it + 1
                    )
                    valid_combos.extend([key + combo for combo in combos])

            return valid_combos


    def __create_dict(
        self, protein: Protein, protein_sequence, keys, depth, step_size, best_options=[], positions=[]
    ):
        """
        Create a dictionary of possible folding sequences and their corresponding scores.

        Parameters:
        - protein (Protein): The protein structure being folded.
        - protein_sequence (str): The protein sequence to be folded.
        - keys (list): The possible folding directions (e.g., ["R", "L", "U", "D"]).
        - depth (int): The depth at which to stop the folding process.
        - step_size (int): The step size used during the folding process.
        - best_options (list): List of previously determined best folding options.

        Returns:
        - dict: A dictionary mapping folding sequences to their scores.
        """
        valid_pos = True


        if len(positions) == 0:
            seq, score_dict, pos = "", {}, [(0, 0, 0)]
        else:
            seq, score_dict, pos = "", {}, [positions[-1]]

        if depth > len(protein_sequence):
            depth = len(protein_sequence)

        valid_combos = self.__valid_combinations(keys, length=depth)

        for steps in valid_combos:
            if len(best_options) != 0:
                for option in best_options:
                    if steps[: depth-step_size] in option and len(steps) > len(option):
                        prt = Protein(protein_sequence[: depth + 1])
                        if len(positions) > 0:
                            prt.get_head().position = positions[-1]

                        dict_ = self.__create_nested_dict(protein, keys, depth + 1)

                        aminoacid_ = prt._head.link

                        for step in steps:
                            dict_ = dict_[step]
                            if isinstance(dict_, dict):
                                aminoacid_.position = dict_["pos"]
                                pos.append(aminoacid_.position)
                                seq += step
                            else:
                                aminoacid_.position = dict_[-1]
                                pos.append(aminoacid_.position)
                                seq += step

                            aminoacid_ = aminoacid_.link

                        current = prt.get_head()
                        while current is not None:
                            if current.position in positions:
                                valid_pos = False
                            prt.add_to_grid(current.position, current)
                            current = current.link

                        if prt.is_valid() and valid_pos:
                            score_dict[seq] = prt.get_score()
                            if len(positions) == 0:
                                seq, pos = "", [(0, 0, 0)]
                            else:
                                seq, pos = "", [positions[-1]]
                        else:
                            current = prt.get_head()
                            while current is not None:
                                prt.remove_from_grid(current.position)
                                current = current.link
                            if len(positions) == 0:
                                seq, pos = "", [(0, 0, 0)]
                            else:
                                seq, pos = "", [positions[-1]]
                            continue
                    else:
                        continue
            else:
                prt = Protein(protein_sequence[: depth + 1])
                if len(positions) > 0:
                            prt.get_head().position = positions[-1]
                dict_ = self.__create_nested_dict(protein, keys, depth + 1)

                aminoacid_ = prt._head.link

                for step in steps:
                    dict_ = dict_[step]
                    if isinstance(dict_, dict):
                        aminoacid_.position = dict_["pos"]
                        pos.append(aminoacid_.position)
                        seq += step
                    else:
                        aminoacid_.position = dict_[-1]
                        pos.append(aminoacid_.position)
                        seq += step

                    aminoacid_ = aminoacid_.link

                current = prt.get_head()
                while current is not None:
                    if current.position in positions:
                        print(current.position, positions)
                        valid_pos = False
                    prt.add_to_grid(current.position, current)
                    current = current.link

                if prt.is_valid() and valid_pos:
                    score_dict[seq] = prt.get_score()
                    if len(positions) == 0:
                        seq, pos = "", [(0, 0, 0)]
                    else:
                        seq, pos = "", [positions[-1]]
                else:
                    current = prt.get_head()
                    while current is not None:
                        prt.remove_from_grid(current.position)
                        current = current.link
                    if len(positions) == 0:
                        seq, pos = "", [(0, 0, 0)]
                    else:
                        seq, pos = "", [positions[-1]]
                    continue

        return score_dict

    def __bfsfold(self, protein: Protein, when_cutting, step, positions) -> Protein:
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
        length_protein = len(protein)
        sequence_protein = protein._sequence
        if self.dimensions == 2:
            types = ["R", "L", "U", "D"]
            when_cutting = 4
        elif self.dimensions == 3:
            types = ["R", "L", "U", "D", "F", "B"]
            when_cutting = 6
        min_keys = []


        step = 1

        for depth in range(when_cutting, length_protein, step):
            create_d = self.__create_dict(
                protein, sequence_protein, types, depth, step, min_keys, positions
            )
            if depth == 4:
                print(create_d)
            print(depth, 1)

            min_key = min(create_d, key=lambda k: create_d[k])
            min_keys = [k for k, v in create_d.items() if v == create_d[min_key]]


        nested_dict = self.__create_nested_dict(protein, types, length_protein)
        aminoacid_ = protein.get_head().link
        
        folding = random.choice(min_keys)

        for direction in folding:

            nested_dict = nested_dict[direction]
            aminoacid_.position = nested_dict["pos"]
            aminoacid_ = aminoacid_.link


        return protein



    def run(self) -> Protein:
        """
        Run the BFS folding algorithm on the specified protein.

        Returns:
        - Protein: The folded protein structure.
        """
        sequences, proteins, positions = self.cutting_seq(self._protein, cut=9), [], []
        for c, seq in enumerate(sequences):
            print(seq)
            if c == 0:
                proteins.append(self.__bfsfold(Protein(seq), self._cut, self._step, positions))
            else:
                for protein in proteins:
                    current = protein.get_head()
                    while current is not None:
                        positions.append(current.position)
                        current = current.link
                    print(positions)
                    proteins.append(self.__bfsfold(Protein(seq), self._cut, self._step, positions))

            positions = []

        for protein in proteins:
            current = protein.get_head()
            positions.append(current.position)
            while current is not None:
                positions.append(current.position)
                current = current.link
        
        prt = Protein(self._sequence)
        print(prt.get_grid())
        current = prt.get_head()

        for pos in positions[1:]:
            current.position = pos
            prt.add_to_grid(current.position, current)
            print(current.position)
            current = current.link

        self._protein, self._sequence = prt, prt._sequence
        start_time = time.time()  # Record the start time


        end_time = time.time()  # Record the end time
        elapsed_time = end_time - start_time
        print(f"Elapsed Time: {elapsed_time} seconds")

        return self._protein

