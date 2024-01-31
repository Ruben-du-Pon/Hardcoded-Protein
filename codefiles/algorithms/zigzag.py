from ..classes.protein import Protein
from typing import List, Tuple


class ZigzagFold:
    """
    Represents a folding algorithm that folds a protein in a zigzag pattern.

    Attributes:
    - movements (List[Tuple[int, int, int]]): The possible movements in the grid.

    Methods:
    - __init__(self, protein: Protein, dimensions: int) -> None: Initializes a
        new instance of the ZigzagFold class.
    - run(self) -> Protein: Runs the zigzag folding algorithm and returns the
        folded protein.
    """  # noqa

    movements: List[Tuple[int, int, int]] = [
        (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, 1, 0)
    ]

    def __init__(self, protein: Protein, dimensions: int) -> None:
        """
        Initializes a new instance of the ZigzagFold class.

        Parameters:
        - protein (Protein): The protein to fold.
        - dimensions (int): The number of dimensions for folding (2 or 3).

        Raises:
        - ValueError: If the dimensions parameter is not 2 or 3.
        """
        if dimensions not in (2, 3):
            raise ValueError(
                "Please enter dimensions as 2 or 3 for 2D or 3D folding.")

        self._protein = protein
        self._dimensions = dimensions

    def run(self) -> Protein:
        """
        Runs the zigzag folding algorithm and returns the folded protein.

        Returns:
        - Protein: The folded protein.

        Raises:
        - ValueError: If a valid folding for the protein cannot be found.
        """
        current = self._protein.get_list()
        self._protein.add_to_grid(current.position, current)

# Initialize the movement index.
        movement_index: int = 0

# Fold the protein until it is valid.
        while not self._protein.is_valid():
            current = current.link

# Add the movement to the position of the previous amino acid.
            current.position = tuple(
                sum(x) for x in zip(current.position,
                                    current.predecessor.position,
                                    ZigzagFold.movements[movement_index])
            )

# Add the amino acid to the grid and increase the movement index,
            # looping back when it reaches the end of the list.
            self._protein.add_to_grid(current.position, current)
            movement_index = (movement_index + 1) % len(ZigzagFold.movements)

        if not self._protein.is_valid():
            raise ValueError(
                f"Couldn't find a valid folding for protein {self._protein}")

        return self._protein
