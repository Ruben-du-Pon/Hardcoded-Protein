from ..classes.protein import Protein
from typing import List, Tuple


class SpiralFold:
    """
    A class that folds a protein into a spiral.

    Attributes
    ----------
    _protein : Protein
        The protein to fold.
    _dimensions : int
        The dimensions of the protein.
    _movements : List[Tuple[int, int, int]]
        The possible movements in the grid.
    _steps : int
        The number of steps to take in the current direction.

    Methods
    -------
    run() -> Protein:
        Runs the spiral folding algorithm.
    """
    movements: List[Tuple[int, int, int]] = [
        (0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0)
    ]

    def __init__(self, protein: Protein, dimensions: int) -> None:
        """
        Initializes a SpiralFold object.

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
        self._steps: int = 1

    def run(self) -> Protein:
        """
        Runs the spiral folding algorithm.

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
        current = self._protein.get_head()
        self._protein.add_to_grid(current.position, current)

        # Initialize a movement index and move to the next amino acid.
        movement_index: int = 0
        current = current.link

        # Keep setting positions until a valid folding is found.
        while not self._protein.is_valid():

            # Set the position of self._steps amino acids in the current direction.
            for _ in range(self._steps):
                current.position = tuple(
                    sum(x) for x in zip(current.predecessor.position,
                                        SpiralFold.movements[movement_index])
                )
                self._protein.add_to_grid(current.position, current)

                # If the end of the protein is reached, stop.
                if current.link is None:
                    break
                current = current.link

            # Change the direction and increase the number of steps.
            movement_index = (movement_index + 1) % len(SpiralFold.movements)
            self._steps += 1 if movement_index % 2 == 0 else 0

        if not self._protein.is_valid():
            raise ValueError(
                f"Couldn't find a valid folding for protein {self._protein}")

        return self._protein
