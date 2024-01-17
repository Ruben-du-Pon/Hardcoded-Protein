from ..classes.protein import Protein
from typing import List, Tuple


class ZigzagFold:
    def __init__(self, protein: Protein, dimensions: int) -> None:
        """
        Initialize the ZigzagFold object.

        Parameters
        ----------
        protein : Protein
            The protein structure to which the zigzag folding algorithm is applied.
        """
        self._protein = protein
        self._dimensions = dimensions
        self._movements: List[Tuple[int, int, int]] = [
            (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, 1, 0)
        ]

    def run(self) -> None:
        """
        Apply a zigzag folding algorithm to the given protein.

        Raises
        ------
        ValueError
            If a valid folding cannot be found for the given protein.
        """
        current = self._protein.get_list()
        self._protein.add_to_grid(current.position, current)
        movement_index: int = 0

        while not self._protein.is_valid():
            current = current.link
            current.position = tuple(
                sum(x) for x in zip(current.position, current.predecessor.position, self._movements[movement_index])
            )
            self._protein.add_to_grid(current.position, current)
            movement_index = (movement_index + 1) % len(self._movements)

        if not self._protein.is_valid():
            raise ValueError(
                f"Couldn't find a valid folding for protein {self._protein}")
