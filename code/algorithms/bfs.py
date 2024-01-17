from ..classes.protein import Protein


class BfsFold:
    def __init__(self, protein: Protein, when_cutting=2, step=1):
        """
        Initialize the BFSFold object.

        Args:
        - protein (Protein): The protein structure to fold.
        - when_cutting (int): The depth at which to start considering alternative folding directions.
        - step (int): The step size for adjusting the depth during the folding process.
        """
        self._protein = protein
        self._cut = when_cutting
        self._step = step

    def __create_nested_dict(
        self, protein: Protein, keys, depth, prev=None, pos=[(0, 0, 0)]
    ):
        """
        Recursively create a nested dictionary representing possible folding paths.

        Args:
        - protein (Protein): The protein structure.
        - keys (list): List of possible directions (e.g., ["R", "L", "U", "D"]).
        - depth (int): The current depth in the recursion.
        - prev (str): The previous direction.
        - t (list): List of tuples representing positions.

        Returns:
        - dict: Nested dictionary representing folding paths.
        """
        aminoacid = protein._head

        if aminoacid.link is None or depth == 1:
            return pos

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
                move = {
                    "R": (1, 0, 0),
                    "L": (-1, 0, 0),
                    "U": (0, 1, 0),
                    "D": (0, -1, 0),
                }
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
        Generate valid combinations of directions for a given length.

        Args:
        - keys (list): List of possible directions (e.g., ["R", "L", "U", "D"]).
        - prev (str): The previous direction.
        - length (int): The length of combinations to generate.
        - it (int): Current iteration.

        Returns:
        - list: List of valid combinations of directions.
        """
        if it == length:
            return [[]]

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
                valid_combos.extend([[key] + combo for combo in combos])

        return valid_combos

    def __create_dict(
        self, protein: Protein, protein_sequence, keys, depth, best_options=[]
    ):
        """
        Create a dictionary of folding paths and their scores.

        Args:
        - protein (Protein): The protein structure.
        - protein_sequence (str): The protein sequence.
        - keys (list): List of possible directions (e.g., ["R", "L", "U", "D"]).
        - depth (int): The maximum depth to explore.
        - best_options (list): List of best folding options.

        Returns:
        - dict: Dictionary mapping folding paths to their scores.
        """
        seq, seq_, score_dict, pos = "", "", {}, [(0, 0, 0)]

        if depth > len(protein_sequence):
            depth = len(protein_sequence)

        for lvl in range(1, depth + 1):
            valid_combos = self.__valid_combinations(keys, length=lvl)
            for steps in valid_combos:
                if len(best_options) != 0:
                    for step in steps:
                        seq_ += step
                    for option in best_options:
                        if seq_[: len(option)] in option and len(seq_) >= len(option):
                            prt = Protein(protein_sequence[: lvl + 1])
                            dict_ = self.__create_nested_dict(protein, keys, lvl + 1)

                            aminoacid_ = prt._head
                            aminoacid_ = aminoacid_.link

                            for step in steps:
                                dict_ = dict_[step]
                                if isinstance(dict_, dict):
                                    aminoacid_.position = dict_["pos"]
                                    pos = pos + [aminoacid_.position]
                                    seq += step
                                else:
                                    aminoacid_.position = dict_[-1]
                                    pos = pos + [aminoacid_.position]
                                    seq += step

                                aminoacid_ = aminoacid_.link

                            current = prt.get_list()
                            prt.add_to_grid(current.position, current)
                            current = current.link
                            while current is not None:
                                prt.add_to_grid(current.position, current)
                                current = current.link

                            if prt.is_valid():
                                score_dict[seq] = prt.get_score()
                                seq, pos = "", [(0, 0, 0)]
                            else:
                                seq, seq_, pos = "", "", [(0, 0, 0)]
                                continue
                        else:
                            seq_ = ""
                            continue
                else:
                    prt = Protein(protein_sequence[: lvl + 1])
                    dict_ = self.__create_nested_dict(protein, keys, lvl + 1)

                    aminoacid_ = prt._head
                    aminoacid_ = aminoacid_.link

                    for step in steps:
                        dict_ = dict_[step]
                        if isinstance(dict_, dict):
                            aminoacid_.position = dict_["pos"]
                            pos = pos + [aminoacid_.position]
                            seq += step
                        else:
                            aminoacid_.position = dict_[-1]
                            pos = pos + [aminoacid_.position]
                            seq += step

                        aminoacid_ = aminoacid_.link

                    current = prt.get_list()
                    prt.add_to_grid(current.position, current)
                    current = current.link
                    while current is not None:
                        prt.add_to_grid(current.position, current)
                        current = current.link

                    if prt.is_valid():
                        score_dict[seq] = prt.get_score()
                        seq, pos = "", [(0, 0, 0)]
                    else:
                        seq, pos = "", [(0, 0, 0)]
                        continue

        return score_dict

    def __bfsfold(self, protein: Protein, when_cutting, step, dimension=2) -> Protein:
        """
        Perform breadth-first folding of a protein structure to optimize its score.

        Args:
        - protein (Protein): The protein structure to fold.
        - when_cutting (int): The depth at which to start considering alternative folding directions.
        - step (int): The step size for adjusting the depth during the folding process.

        Returns:
        - Protein: The folded protein structure.
        """
        length_protein = len(protein)
        sequence_protein = protein._sequence
        types = ["R", "L", "U", "D"]
        min_keys = []

        when_cutting = length_protein // 2

        for depth in range(when_cutting, length_protein, step):
            create_d = self.__create_dict(
                protein, sequence_protein, types, depth, min_keys
            )
            min_key = min(create_d, key=lambda k: create_d[k])
            min_keys = [k for k, v in create_d.items() if v == create_d[min_key]]

        nested_dict = self.__create_nested_dict(protein, types, length_protein)
        aminoacid_ = protein.get_list()

        for direction in min_keys[-1]:
            nested_dict = nested_dict[direction]
            aminoacid_ = aminoacid_.link
            if isinstance(nested_dict, dict):
                aminoacid_.position = nested_dict["pos"]
            else:
                aminoacid_.position = nested_dict[-1]

        aminoacid_ = protein._head

        while aminoacid_ is not None:
            protein.add_to_grid(aminoacid_.position, aminoacid_)
            aminoacid_ = aminoacid_.link

        return protein

    def run(self) -> Protein:
        """
        Run the BFS folding algorithm on a given Protein instance.

        Returns:
        - Protein: The folded protein structure.
        """
        return self.__bfsfold(self._protein, self._cut, self._step)
