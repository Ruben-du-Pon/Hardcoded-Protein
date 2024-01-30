from ..classes.protein import Protein
from typing import List, Tuple


class ZigzagFold:
    """
    A class that folds a protein into a zigzag.

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
    run() -> None
        Runs the zigzag folding algorithm.
    """
    # Define the possible movements in the grid.
    movements: List[Tuple[int, int, int]] = [
        (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, 1, 0)
    ]

    def __init__(self, protein: Protein, dimensions: int) -> None:
        """
        Initializes a ZigzagFold object.

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

    def run(self) -> Protein:
        """
        Runs the zigzag folding algorithm.

        Returns
        -------
        Protein
            The folded protein.

        Raises
        ------
        ValueError
            If a valid folding cannot be found for the given protein.
        """
        # Add the first amino acid to the grid.
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
