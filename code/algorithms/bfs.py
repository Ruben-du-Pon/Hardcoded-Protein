from ..classes.protein import Protein


class BfsFold:
    def __init__(self, protein: Protein, when_cutting=2, step=1):
        self._protein = protein
        self._cut = when_cutting
        self._step = step

    def __create_nested_dict(
        self, protein: Protein, keys, depth, prev=None, pos=[(0, 0, 0)]
    ):
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

    def __valid_combinations(self, keys, prev=None, length=2, it=0):
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
        seq, seq_, score_dict, pos = "", "", {}, [(0, 0, 0)]

        if depth > len(protein_sequence):
            depth = len(protein_sequence)

        valid_combos = self.__valid_combinations(keys, length=depth)

        for steps in valid_combos:
            if len(best_options) != 0:
                for step in steps:
                    seq_ += step
                for option in best_options:
                    if seq_[: len(option)] in option and len(seq_) > len(option):
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
                            seq, seq_, pos = "", "", [(0, 0, 0)]
                            continue
                    else:
                        seq_ = ""
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
                    seq, pos = "", [(0, 0, 0)]
                    continue

        return score_dict

    def __bfsfold(self, protein: Protein, when_cutting, step, dimension=2) -> Protein:
        length_protein = len(protein)
        sequence_protein = protein._sequence
        types = ["R", "L", "U", "D"]
        min_keys = []

        when_cutting = 6
        step = 2

        for depth in range(when_cutting, length_protein, step):
            create_d = self.__create_dict(
                protein, sequence_protein, types, depth, min_keys
            )

            min_key = min(create_d, key=lambda k: create_d[k])
            min_keys = [k for k, v in create_d.items() if v == create_d[min_key]]
            max_length = max(len(key) for key in min_keys)

            if max_length <= 11:
                min_keys = [key for key in min_keys if len(key) == max_length]
            else:
                min_keys = [min_keys[-1]]

        nested_dict = self.__create_nested_dict(protein, types, length_protein)
        aminoacid_ = protein.get_head().link

        for direction in min_keys[-1]:
            nested_dict = nested_dict[direction]
            aminoacid_.position = nested_dict["pos"]
            aminoacid_ = aminoacid_.link

        aminoacid_ = protein.get_head().link

        while aminoacid_ is not None:
            protein.add_to_grid(aminoacid_.position, aminoacid_)
            aminoacid_ = aminoacid_.link

        return protein

    def run(self) -> Protein:
        return self.__bfsfold(self._protein, self._cut, self._step)
