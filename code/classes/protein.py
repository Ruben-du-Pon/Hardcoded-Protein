from code.classes.aminoacid import Aminoacid
from operator import add, sub
import csv
import random


class Protein:
    def __init__(self, sequence: str) -> None:
        self._sequence: str = sequence
        self._length: int = len(sequence)
        # Store positions in the instance
        self._positions: list[tuple[int]
                              ] = self.generate_random_protein_positions()
        self._head: Aminoacid = self.create_double_linked_list(self)
        self._grid: dict[tuple, Aminoacid] = {}
        self._score = 0

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

        head = Aminoacid(type=self._sequence[0], position=self._positions[0])
        current = head

        for idx, type in enumerate(self._sequence[1:], start=1):
            new_aminoacid = Aminoacid(
                type=type, position=self._positions[idx], predecessor=current)
            current.link = new_aminoacid
            current = new_aminoacid

        return head

    def get_score(self):
        return self._score

    def get_folding(self):
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
                """
                x, y, z = current.position
                next_x, next_y, next_z = next_amino.position

                # Calculate the direction based on position change
                if x != next_x:
                    fold = 1 if next_x > x else -1
                elif y != next_y:
                    fold = 2 if next_y > y else -2
                elif z != next_z:
                    fold = 3 if next_z > z else -3
                else:
                    fold = 0  # No change in position, might need to handle differently"""

            folding.append({'amino': current.get_type(), 'fold': fold})
            current = current.link

        folding.append({'amino': 'score', 'fold': self._score})
        return folding

    def create_csv(self, index: int = 0) -> None:
        """
        Creates a csv file that displays a specific folding of a protein
        post: creates output.csv if it doesn't exist, empties it if it does,
        then fills it with the folding data.

        pre: index is an int, the get_folding method outputs a list of dicts with
        keys amino and fold, and resp values P, H or C and 1, -1, 2, -2, 3 or -3.
        post: creates output[index].csv with a header amino, score; a footer score, <score>
        and a body with P, H or C followed by direction 1, -1, 2, -2, 3 or -3.
        """

        filename = "data/output/output" + str(index) + ".csv"
        with open(filename, 'w', newline='') as file:
            print(filename)
            header = ["amino", "fold"]
            writer = csv.DictWriter(file, fieldnames=header)
            folding = self.get_folding()

            writer.writeheader()
            writer.writerows(folding)

    def generate_random_protein_positions(self):
        directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
        positions = [(0, 0, 0)]  # Starting position

        for _ in range(1, self._length):
            last_position = positions[-1]
            movement = random.choice(directions)
            next_position = tuple(map(add, last_position, movement))
            positions.append(next_position)

        return positions

    def get_grid_2D(self):
        # Gather sequence of amino acids directly in the method
        sequence = []
        current = self._head
        while current:
            sequence.append(current.get_type())
            current = current._link

        # Determine the grid size
        min_x = min(x for x, _, _ in self._positions)
        max_x = max(x for x, _, _ in self._positions)
        min_y = min(y for _, y, _ in self._positions)
        max_y = max(y for _, y, _ in self._positions)
        grid_size_x = (max_x - min_x) * 2 + 3
        grid_size_y = (max_y - min_y) * 2 + 3

        # Create an empty grid
        grid = [[' ' for _ in range(grid_size_x)] for _ in range(grid_size_y)]

        # Adjust positions for the grid
        adjusted_positions = [(2 * (x - min_x) + 1, 2 * (y - min_y) + 1)
                              for x, y, _ in self._positions]

        # Place amino acids and connections
        for (x, y), amino_acid in zip(adjusted_positions, sequence):
            grid[y][x] = '' + amino_acid + ''
        for (x1, y1), (x2, y2) in zip(adjusted_positions, adjusted_positions[1:]):
            if x1 == x2:
                for y in range(min(y1, y2) + 1, max(y1, y2)):
                    grid[y][x1] = '|'
            elif y1 == y2:
                for x in range(min(x1, x2) + 1, max(x1, x2)):
                    grid[y1][x] = '-'

        # Convert the grid to a string representation
        grid_str = '\n'.join([''.join(row) for row in grid])
        return grid_str


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
"""
