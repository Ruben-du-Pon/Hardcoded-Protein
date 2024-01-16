from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid

"""
NOTE: Still busy working on the algotihm.
NOTE: ONLY WORKS FOR 2D FOR NOW.
NOTE: A bit slow for now.
"""

def create_nested_dict(protein: Protein, keys, depth, prev=None, t=[(0, 0, 0)]):
    """
    Recursively create a nested dictionary representing possible protein folding steps.

    Args:
    - protein (Protein): The protein structure.
    - keys (list): The possible directions to fold the protein (e.g., ["R", "L", "U", "D"]).
    - depth (int): The depth of recursion.
    - prev (str): The previous direction of folding.
    - t (list): A list of tuples representing the current position in the folding process.

    Returns:
    - dict: A nested dictionary representing possible protein folding steps.
    """
    aminoacid = protein._head

    if aminoacid.link is None or depth == 1:
        return t

    result_dict = {"pos": t[-1]}
    for key in keys:
        if (key == "R" and prev == "L") or \
                (key == "L" and prev == "R") or \
                (key == "U" and prev == "D") or \
                (key == "D" and prev == "U"):
            continue
        else:
            move = {"R": (1, 0, 0), "L": (-1, 0, 0), "U": (0, 1, 0), "D": (0, -1, 0)}
            result_dict[key] = create_nested_dict(protein, keys, depth - 1, key, t + [tuple(x + y for x, y in zip(t[-1], move[key]))])

    return result_dict


def valid_combinations(keys, prev=None, length=2, it=0):
    """
    Generate valid combinations of folding directions up to a specified length.

    Args:
    - keys (list): The possible directions to fold the protein (e.g., ["R", "L", "U", "D"]).
    - prev (str): The previous direction of folding.
    - length (int): The desired length of folding directions.
    - it (int): Current iteration level.

    Returns:
    - list: List of valid combinations of folding directions.
    """
    if it == length:
        return [[]]

    valid_combos = []

    for key in keys:
        if (key == "R" and prev == "L") or \
           (key == "L" and prev == "R") or \
           (key == "U" and prev == "D") or \
           (key == "D" and prev == "U"):
            continue
        else:
            combos = valid_combinations(keys, prev=key, length=length, it=it+1)
            valid_combos.extend([[key] + combo for combo in combos])

    return valid_combos


def create_dict(protein:Protein, protein_sequence, keys, depth, best_options=[]):
    """
    Create a dictionary of protein sequences and their scores based on valid folding directions.

    Args:
    - protein (Protein): The protein structure.
    - protein_sequence (str): The sequence of amino acids in the protein.
    - keys (list): The possible directions to fold the protein (e.g., ["R", "L", "U", "D"]).
    - depth (int): The maximum depth of folding directions to consider.
    - best_options (list): List of optimal folding sequences from previous iterations.

    Returns:
    - dict: A dictionary of protein sequences and their corresponding scores.
    """
    seq, seq_ , score_dict, pos = "", "", {}, [(0,0,0)]

    if depth > len(protein_sequence):
        depth = len(protein_sequence)

    for lvl in range(1, depth+1):
        valid_combos = valid_combinations(keys, length=lvl)
        for steps in valid_combos:
            if len(best_options) != 0:
                for step in steps:
                    seq_ += step
                for option in best_options:
                    if seq_[:len(option)] in option and len(seq_) >= len(option):
                        
                        prt = Protein(protein_sequence[:lvl+1])
                        dict_ = create_nested_dict(protein, keys, lvl+1)

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
                        
                        if len(pos) == len(set(pos)):
                            score_dict[seq] = prt.calculate_score()
                            seq, seq_, pos = "", "", [(0,0,0)]
                        else:
                            seq, seq_, pos = "", "", [(0,0,0)]
                            continue
                    else:
                        seq_ = ""
                        continue
            else:
                prt = Protein(protein_sequence[:lvl+1])
                dict_ = create_nested_dict(protein, keys, lvl+1)

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

                if len(pos) == len(set(pos)):
                    score_dict[seq] = prt.calculate_score()
                    seq, pos = "", [(0,0,0)]
                else:
                    seq, pos = "", [(0,0,0)]
                    continue


    return score_dict


def breadth_fold(protein: Protein, when_cutting=4, step=2) -> Protein:
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


    for depth in range(when_cutting, length_protein, step):
        create_d = create_dict(protein, sequence_protein, types, depth, min_keys)
        min_key = min(create_d, key=lambda k: create_d[k])
        min_keys = [k for k, v in create_d.items() if v == create_d[min_key]]
    
    nested_dict = create_nested_dict(protein, types, length_protein)

    # protein_ = Protein(protein._sequence)
    aminoacid_ = protein.get_list()

    for direction in min_key:
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

