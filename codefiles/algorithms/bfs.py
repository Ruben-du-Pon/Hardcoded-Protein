from ..classes.protein import Protein
from ..classes.protein import Aminoacid
import random
from typing import List, Set, Tuple
import time
import numpy as np


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

    def cutting_seq(self, protein: Protein, cut=8) -> Set[str]:
        length, seq, proteins = len(protein), protein._sequence, []

        for x in range(0, length, cut):
            proteins.append(Protein(seq[x:x+cut]))

        return proteins

    def __valid_combinations(self, keys, prev=None, length=2, it=0, prev_valid=set()) -> List[str]:
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
                return {''}

            valid_combos = set()

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
                    valid_combos.update({key + combo for combo in combos})

            return valid_combos

        elif self.dimensions == 3:
            if it == length:
                return {''}

            valid_combos = set()

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
                    valid_combos.update({key + combo for combo in combos})

            return valid_combos

    def __add_combinations(self, prev_valid):
        new_foldings = set()

        if self.dimensions == 2:
            types = {"R", "L", "U", "D"}
        elif self.dimensions == 3:
            types = {"R", "L", "U", "D", "F", "B"}

        for prev in prev_valid:
            last = prev[-1]
            if self.dimensions == 2 or self.dimensions == 3:
                if last == "R":
                    types.remove("L")
                elif last == "L":
                    types.remove("R")
                elif last == "U":
                    types.remove("D")
                elif last == "D":
                    types.remove("U")
            if self.dimensions == 3:
                if last == "F":
                    types.remove("B")
                elif last == "B":
                    types.remove("F")

            for direction in types:
                new_foldings.add(prev+direction)

            if self.dimensions == 2:
                types = {"R", "L", "U", "D"}
            elif self.dimensions == 3:
                types = {"R", "L", "U", "D", "F", "B"}

        return new_foldings

    def __create_nested_dict(
        self, protein: Protein, keys, depth, prev=None, pos=[(0, 0, 0)]
    ) -> dict:
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
                    move = {"R": (1, 0, 0), "L": (-1, 0, 0),
                            "U": (0, 1, 0), "D": (0, -1, 0)}
                    result_dict[key] = self.__create_nested_dict(
                        protein,
                        keys,
                        depth - 1,
                        key,
                        pos +
                        [tuple(x + y for x, y in zip(pos[-1], move[key]))],
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
                    move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0),
                            "D": (0, -1, 0), "F": (0, 0, 1), "B": (0, 0, -1)}
                    result_dict[key] = self.__create_nested_dict(
                        protein,
                        keys,
                        depth - 1,
                        key,
                        pos +
                        [tuple(x + y for x, y in zip(pos[-1], move[key]))],
                    )

            return result_dict

    def __create_dict(
        self, protein: Protein, protein_sequence, keys, depth, step_size, best_options=set(), posit=None
    ) -> dict:
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
        seq, score_dict, pos = "", {}, [(0, 0, 0)]

        if depth > len(protein_sequence):
            depth = len(protein_sequence)

            if best_options != set():
                valid_combos = self.__add_combinations(best_options)
            else:
                valid_combos = self.__valid_combinations(keys, length=depth)

            for steps in valid_combos:
                if len(best_options) != 0:
                    for option in best_options:
                        # print(9, best_options)
                        # print(10, valid_combos)
                        if option in steps[: depth - step_size+1] and len(steps) > len(option):
                            prt = Protein(protein_sequence[: depth + 1])

                            current = prt.get_head()
                            for i in posit:
                                current.position = i
                                if i != (0, 0, 0):
                                    pos.append(i)
                                current = current.link

                            dict_ = self.__create_nested_dict(
                                protein, keys, 2, pos=[posit[-1]])

                            # print(dict_)

                            for step in steps[-1]:
                                dict_ = dict_[step]
                                if isinstance(dict_, dict):
                                    current.position = dict_["pos"]
                                    pos.append(current.position)
                                    seq += step
                                else:
                                    current.position = dict_[-1]
                                    pos.append(current.position)
                                    seq += step

                            current = prt.get_head()
                            while current is not None:
                                # print(current.position)
                                prt.add_to_grid(current.position, current)
                                current = current.link

                            if prt.is_valid():
                                score_dict[option+seq] = prt.get_score()
                                seq, pos = "", [(0, 0, 0)]
                            else:
                                aminoacid_.position = dict_[-1]
                                pos.append(aminoacid_.position)
                                seq += step

                            aminoacid_ = aminoacid_.link

                        current = prt.get_head()
                        while current is not None:
                            prt.add_to_grid(current.position, current)
                            current = current.link

                        if prt.is_valid():
                            score_dict[seq] = prt.get_score()
                            seq, pos = "", [(0, 0, 0)]
                        else:
                            current = prt.get_head()
                            while current is not None:
                                prt.remove_from_grid(current.position)
                                current = current.link
                            seq, pos = "", [(0, 0, 0)]
                            continue
                    else:
                        continue
            else:
                prt = Protein(protein_sequence[: depth + 1])
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
                    prt.add_to_grid(current.position, current)
                    current = current.link

                if prt.is_valid():
                    score_dict[seq] = prt.get_score()
                    seq, pos = "", [(0, 0, 0)]
                else:
                    current = prt.get_head()
                    while current is not None:
                        prt.remove_from_grid(current.position)
                        current = current.link
                    seq, pos = "", [(0, 0, 0)]
                    continue

        return score_dict

    def __bfsfold(self, protein: Protein, when_cutting, step) -> Protein:
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
        posit = [(0, 0, 0)]
        sequence_protein = protein._sequence
        if self.dimensions == 2:
            types = {"R", "L", "U", "D"}
            move = {"R": (1, 0, 0), "L": (-1, 0, 0),
                    "U": (0, 1, 0), "D": (0, -1, 0)}
        elif self.dimensions == 3:
            types = {"R", "L", "U", "D", "F", "B"}
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0),
                    "D": (0, -1, 0), "F": (0, 0, 1), "B": (0, 0, -1)}
        min_keys = set()
        
        if "C" in self._sequence:
            when_cutting = 6
        else:
            when_cutting = 4
        
        step = 1

        for depth in range(when_cutting, length_protein, step):
            # print(8, min_keys)
            create_d = self.__create_dict(
                protein, sequence_protein, types, depth, step, min_keys, posit
            )
            posit = [(0, 0, 0)]

            min_key = min(create_d, key=lambda k: create_d[k])
            min_keys = {k for k, v in create_d.items() if v ==
                        create_d[min_key]}
            # print(2, min_keys)
            if len(min_keys) >= 2:
                # Randomly select 1 of the set
                min_keys = random.sample(min_keys, 1)
            else:
                min_keys = list(min_keys)
            for key in min_keys[0]:
                posit.append(tuple(np.array(posit[-1]) + np.array(move[key])))

        nested_dict = self.__create_nested_dict(protein, types, length_protein)
        aminoacid_ = protein.get_head().link

        folding = random.choice(list(min_keys))

        for direction in folding:
            nested_dict = nested_dict[direction]
            aminoacid_.position = nested_dict["pos"]
            aminoacid_ = aminoacid_.link

        aminoacid_ = protein.get_head()

        while aminoacid_ is not None:
            protein.add_to_grid(aminoacid_.position, aminoacid_)
            aminoacid_ = aminoacid_.link

        return protein

    def cut_protein_seq(self, protein: Protein, dimensions=2) -> List[Protein]:
        seq, seq_len, lst_proteins = protein._sequence, len(protein), []

        if dimensions == 2:
            for seq_l in range(0, seq_len, 10):
                # print(seq[seq_l:seq_l+10])
                lst_proteins.append(Protein(seq[seq_l:seq_l+10]))
        elif dimensions == 3:
            for seq_l in range(0, seq_len, 8):
                # print(seq[seq_l:seq_l+8])
                lst_proteins.append(Protein(seq[seq_l:seq_l+8]))

        return lst_proteins

    """
    NOTE: This function is for 'Simulated Annealing'
    """

    def get_possible_foldings(self, protein: Protein,
                              first_coordinate: Tuple[int, int, int],
                              last_coordinate: Tuple[int, int, int]) -> \
            List[List[Aminoacid]]:
        valid_foldings, proteins, aminoacids, protein_aminoacids, length = [
        ], [], [], [], len(protein)
        move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0),
                "D": (0, -1, 0), "F": (0, 0, 1), "B": (0, 0, -1)}

        if self.dimensions == 2:
            types = {"R", "L", "U", "D"}
        elif self.dimensions == 3:
            types = {"R", "L", "U", "D", "F", "B"}

        possible_foldings = self.__valid_combinations(
            keys=types, prev=None, length=length-1)
        dict_fold = self.__create_nested_dict(
            protein, types, length, prev=None, pos=[first_coordinate])

        for folding in possible_foldings:
            dict_ = dict_fold
            for step in folding:
                dict_ = dict_[step]

            if dict_['pos'] == last_coordinate:
                valid_foldings.append(folding)

        for folding in valid_foldings:
            prt = Protein(protein._sequence)
            current = prt._head
            current.position = first_coordinate

            for direction in folding:
                current = current.link
                current.position = tuple(
                    np.array(current.predecessor.position) + np.array(move[direction]))

            proteins.append(prt)

        for prt in proteins:
            current = prt._head
            aminoacids.append(current)
            current = current.link
            while current is not None:
                aminoacids.append(current)
                current = current.link

            protein_aminoacids.append(aminoacids)
            aminoacids = []

        return protein_aminoacids

    def run(self) -> Protein:
        """
        Run the BFS folding algorithm on the specified protein.

        Returns:
        - Protein: The folded protein structure.
        """
        start_time = time.time()

        result = self.__bfsfold(self._protein, self._cut, self._step)

        end_time = time.time()  # Record the end time
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time} seconds")

        return result
