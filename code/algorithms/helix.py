from ..classes.protein import Protein
from typing import List, Tuple


class HelixFold:
    """
    A class that folds a protein into a helix.

    Attributes
    ----------
    _protein : Protein
        The protein to fold.
    _dimensions : int
        The dimensions of the protein.
    _movements : List[Tuple[int, int, int]]
        The possible movements in the grid.

    Methods
    -------
    run()
        Runs the helix folding algorithm.
    """

    def __init__(self, protein: Protein, dimensions: int) -> None:
        """
        Initializes a HelixFold object.

        Parameters
        ----------
        protein : Protein
            The protein to fold.
        dimensions : int
            The dimensions of the protein.

        Raises
        ------
        ValueError
            If the dimensions are not 2 or 3.
        """
        if dimensions not in (2, 3):
            raise ValueError(
                "Please enter dimensions as 2 or 3 for 2D or 3D folding.")

        self._protein = protein
        self._dimensions = dimensions
        self._movements: List[Tuple[int, int, int]] = [
            (0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0), (0, 0, 1)
        ]

    def run(self) -> None:
        """
        Runs the helix folding algorithm.

        Raises
        ------
        ValueError
            If a valid folding cannot be found for the given protein.
        """
        current = self._protein.get_list()
        self._protein.add_to_grid(current.position, current)

        movement_index: int = 0
        counter: int = 1

        current = current.link

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
