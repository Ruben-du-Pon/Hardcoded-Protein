from ..classes.protein import Protein
from typing import List, Tuple


class HelixFold:
    def __init__(self, protein: Protein, dimensions: int):
        """
        Initialize the HelixFold object.

        Parameters
        ----------
        protein : Protein
            The protein structure to which the helical folding algorithm is applied.
        dimensions : int
            Represents the dimensions of folding (2D or 3D).
        """

        self._protein = protein
        self._movements: List[Tuple[int, int, int]] = [
            (0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0), (0, 0, 1)
        ]

    def run(self) -> None:
        """
        Apply a helical folding algorithm to the given protein.

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
        counter: int = 1

        while not self._protein.is_valid():
            movement = self._movements[movement_index]
            if counter == 4:
                movement = (0, 0, 1)
                movement_index -= 1
                counter = 0

            current.position = tuple(
                sum(x) for x in zip(current.predecessor.position, movement)
            )
            self._protein.add_to_grid(current.position, current)

            if current.link is None:
                break

            current = current.link
            movement_index = (movement_index + 1) % len(self._movements)
            counter += 1

        if not self._protein.is_valid():
            raise ValueError(
                f"Couldn't find a valid folding for protein {self._protein}"
            )
