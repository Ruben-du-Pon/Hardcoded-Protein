from .aminoacid import Aminoacid
from operator import add, sub
import csv
import random
from typing import Tuple, Dict, List


class Protein:
    def __init__(self, sequence: str) -> None:
        self._sequence: str = sequence
        self._length: int = len(sequence)
        self._grid: Dict[Tuple[int, int], Aminoacid] = {}
        # Store positions in the instance
        # self._positions: List[Tuple[int]
        #                      ] = self.generate_random_protein_positions(backtracking=True)
        self._head: Aminoacid = self.create_double_linked_list(self)
        self._score: int = 0

    @staticmethod
    def create_double_linked_list(self):
        """
        This particular method can't be called for a created object.
        It's used only once, when creating the object itself.
        After the object is created, this particular function can't
        be called on top of it.
        """
        if not self._sequence:
            return None

        head = Aminoacid(type=self._sequence[0])
        current = head

        for idx, type in enumerate(self._sequence[1:], start=1):
            new_aminoacid = Aminoacid(
                type=type, predecessor=current)
            current.link = new_aminoacid
            current = new_aminoacid

        return head

    def get_score(self):
        return self._score

    def get_folding(self) -> List[Dict[str, int]]:
        folding = []
        current = self._head

        if current is None or current.link is None:
            return folding  # Return empty list in this case

        while current is not None:
            # If current is the last amino acid, there's no next position
            if current.link is None:
                fold = 0
            else:
                next_amino = current.link
                difference = tuple(
                    map(sub, next_amino.position, current.position))
                fold = difference[0] + 2 * difference[1] + 3 * difference[2]

            folding.append({'amino': current.get_type(), 'fold': fold})
            current = current.link

        folding.append({'amino': 'score', 'fold': self._score})
        return folding

    def create_csv(self, index: int = 0) -> None:
        """
        Creates a csv file that displays a specific folding of a protein
        post: creates output.csv if it doesn't exist, empties it if it does,
        then fills it with the folding data.

        pre: index is an int, the get_folding method outputs a list of dicts 
        with keys amino and fold, and resp values P, H or C and 1, -1, 2, -2, 3 or -3.
        post: creates output[index].csv with a header amino, score; a footer 
        score, <score> and a body with P, H or C followed by direction 1, -1, 2, -2, 3 or -3.
        """  # noqa

        filename = "data/output/output" + str(index) + ".csv"
        with open(filename, 'w', newline='') as file:
            print(f"{filename} created.")
            header = ["amino", "fold"]
            writer = csv.DictWriter(file, fieldnames=header)
            folding = self.get_folding()

            writer.writeheader()
            writer.writerows(folding)

    def is_valid(self) -> bool:
        """
        Checks if every aminoacid has a different position.

        post: returns True if every aminoacid has a different position,
        False otherwise.
        """
        return self._length == len(self._grid)

    def is_valid_fold(self, position: Tuple[int]) -> bool:
        """
        Checks if a fold is valid.

        post: returns True if the chosen position is not yet a key in self._grid,
        False otherwise.
        """  # noqa
        return position not in self._grid

    def get_list(self) -> Aminoacid:
        return self._head

    def add_to_grid(self, position: Tuple[int], acid: Aminoacid) -> None:
        self._grid[position] = acid

    def __str__(self) -> str:
        return self._sequence


""" Example usage:
protein_sequence, new_seq = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH", []
print(f"\nOriginal Protein-string: {protein_sequence}\n\n")
protein = Protein(protein_sequence)
current_node = protein._head

print("From left to right: \n")
while current_node:
    print(current_node._type, end=" ")
    if current_node._link == None:
        print("\n\n")
        break
    else:
        current_node = current_node._link

print("From right to left: \n")
while current_node:
    print(current_node._type, end=" ")
    new_seq.append(current_node._type)
    if len(protein_sequence) == len(new_seq):
        print("\n\n")
        break
    else:
        current_node = current_node._predecessor

sample_protein = Protein("HHPHHHPHPHHHPH")
print(sample_protein.get_folding())
print(protein.get_grid_2D())
"""  # noqa
