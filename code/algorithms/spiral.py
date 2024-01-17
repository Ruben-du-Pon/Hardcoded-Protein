from ..classes.protein import Protein
from typing import List, Tuple


class SpiralFold:
    def __init__(self, protein: Protein, dimensions: int) -> None:
        """
        Initialize the SpiralFold object.

        Parameters
        ----------
        protein : Protein
            The protein structure to which the spiral folding algorithm is applied.
        dimensions : int
            Represents the dimensions of folding (2 or 3).
        """
        self._protein = protein
        self._dimensions = dimensions
        self._movements: List[Tuple[int, int, int]] = [
            (0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0)
        ]
        self._steps: int = 1

    def run(self) -> None:
        """
        Apply a spiral folding algorithm to the given protein.

        Raises
        ------
        ValueError
            If a valid folding cannot be found for the given protein.
        """
        current = self._protein.get_list().link
        self._protein.add_to_grid(
            current.predecessor.position, current.predecessor)
        self._protein.add_to_grid(current.position, current)
        movement_index: int = 0

        while not self._protein.is_valid():

            for _ in range(self._steps):
                current.position = tuple(
                    sum(x) for x in zip(current.predecessor.position, self._movements[movement_index])
                )
                self._protein.add_to_grid(current.position, current)
                if current.link is None:
                    break
                current = current.link

            movement_index = (movement_index + 1) % len(self._movements)
            self._steps += 1 if movement_index % 2 == 0 else 0

        if not self._protein.is_valid():
            raise ValueError(
                f"Couldn't find a valid folding for protein {self._protein}")
