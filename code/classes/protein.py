from aminoacid import Aminoacid
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
        self._positions: List[Tuple[int]] = self.generate_random_protein_positions()
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
            print(filename)
            header = ["amino", "fold"]
            writer = csv.DictWriter(file, fieldnames=header)
            folding = self.get_folding()

            writer.writeheader()
            writer.writerows(folding)

    def generate_random_protein_positions(self, backtracking=False):
        """
        Generates a random, valid sequence of positions for the amino acids in the protein.

        It uses backtracking to avoid overlapping positions, ensuring a valid 3D structure. 
        Raises ValueError if a valid structure cannot be generated.

        Returns a list of 3D coordinates for each amino acid.
        """  # noqa

        directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
        positions = [(0, 0, 0)]  # Starting position

        if backtracking:
            self._grid[positions[0]] = self._sequence[0]

            if not self._place_next_amino_acid(positions, 1, directions):
                raise ValueError(
                    "Cannot generate a valid protein structure with the given sequence.")  # noqa

            return positions

        for _ in range(1, self._length):
            last_position = positions[-1]
            movement = random.choice(directions)
            next_position = (last_position[0] + movement[0],
                             last_position[1] + movement[1],
                             last_position[2] + movement[2])
            positions.append(next_position)

        return positions

    def _place_next_amino_acid(self, positions, index, directions):
        """
        Recursive method to place an amino acid in a valid position.

        Args:
        positions: List of current amino acid positions.
        index: Index of the amino acid to be placed.
        directions: Possible directions to move from the current position.

        Returns a bool: True if a valid position is found, False otherwise (backtracks if False).
        """  # noqa

        if index >= self._length:
            return True

        last_position = positions[-1]
        random.shuffle(directions)  # Randomize the directions to try

        for movement in directions:
            next_position = tuple(map(add, last_position, movement))

            if self.is_valid_fold(next_position):
                positions.append(next_position)
                self._grid[next_position] = self._sequence[index]

                if self._place_next_amino_acid(positions, index + 1, directions):  # noqa
                    return True

                # Backtrack
                positions.pop()
                del self._grid[next_position]

        return False

    def get_grid_2D(self) -> str:
        """
        Generates a 2D grid representation of the protein structure.

        Returns a string representation of the 2D grid, showing amino acids
        and their connections.
        """

        sequence = []
        current = self._head
        while current:
            sequence.append(current.get_type())
            current = current.link

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
        for (x1, y1), (x2, y2) in zip(adjusted_positions, adjusted_positions[1:]):  # noqa
            if x1 == x2:
                for y in range(min(y1, y2) + 1, max(y1, y2)):
                    grid[y][x1] = '|'
            elif y1 == y2:
                for x in range(min(x1, x2) + 1, max(x1, x2)):
                    grid[y1][x] = '-'

        # Convert the grid to a string representation
        grid_str = '\n'.join([''.join(row) for row in grid])
        return grid_str

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
