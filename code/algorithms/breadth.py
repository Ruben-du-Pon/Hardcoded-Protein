from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid


def create_nested_dict(protein: Protein, keys, depth, t=[(0, 0, 0)]):
    aminoacid = protein._head

    if aminoacid.link is None or depth == 1:
        return t

    result_dict = {}
    for key in keys:
        if key == "R":
            result_dict["pos"] = t[-1]
            result_dict[key] = create_nested_dict(protein, keys, depth - 1, t+[tuple(x + y for x, y in zip(t[-1], (1,0,0)))])
        elif key == "L":
            # print(t, key)
            result_dict["pos"] = t[-1]
            result_dict[key] = create_nested_dict(protein, keys, depth - 1, t+[tuple(x + y for x, y in zip(t[-1], (-1,0,0)))])
        elif key == "U":
            result_dict["pos"] = t[-1]
            result_dict[key] = create_nested_dict(protein, keys, depth - 1, t+[tuple(x + y for x, y in zip(t[-1], (0,1,0)))])
        else:
            result_dict["pos"] = t[-1]
            result_dict[key] = create_nested_dict(protein, keys, depth - 1, t+[tuple(x + y for x, y in zip(t[-1], (0,-1,0)))])

    return result_dict


def breadth_search(protein: Protein, aminoacid: Aminoacid, dictionary, types, seq: str = ""):
    result_dict = {}

    if "pos" not in dictionary:
        aminoacid.position = dictionary[-1]
        return result_dict
    else:
        aminoacid.position = dictionary["pos"]

    score = protein.calculate_score()
    result_dict[seq] = score
    print(seq, protein._head.position, protein._head.link.position, protein._head.link.link.position, protein._head.link.link.position)
    for type in types:
        if type in dictionary.keys():
            result_dict.update(breadth_search(protein, aminoacid.link, dictionary[type], types, seq + type))

    return result_dict

# types = ["R", "L", "U", "D"]
# protein = Protein("HHHHHP")
# depth = 3
# # begin_pos = (0, 0)
# result = create_nested_dict(protein, types, depth)
# print(result)
# result_ = breadth_search(protein, protein._head, result, types)
# print(result_)
# score = 0
# for i in result_:
#     if len(i) == 3:
#         score += 1

# print(score)


# prt = Protein("HHHHP")

# prt._head.position = (0, 0, 0)
# prt._head.link.position = (1, 0, 0)
# prt._head.link.link.position = (1, 1, 0)
# prt._head.link.link.link.position = (0, 1, 0)

# print(prt.calculate_score())