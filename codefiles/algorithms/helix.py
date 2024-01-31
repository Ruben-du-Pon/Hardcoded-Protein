from ..classes.protein import Protein
from typing import List, Tuple


class HelixFold:
    """
    Represents a folding algorithm for proteins in a helix pattern.

    Attributes:
    - movements (List[Tuple[int, int, int]]): List of possible movements in the grid.

    Methods:
    - __init__(self, protein: Protein, dimensions: int) -> None: Initializes a
        new instance of the HelixFold class.
    - run(self) -> Protein: Runs the helix folding algorithm and returns the
        folded protein.

    Raises:
    - ValueError: If the dimensions parameter is not 2 or 3.
    """  # noqa

    movements: List[Tuple[int, int, int]] = [
        (0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0), (0, 0, 1)
    ]

    def __init__(self, protein: Protein, dimensions: int) -> None:
        """
        Initializes a new instance of the HelixFold class.

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
        Runs the helix folding algorithm and returns the folded protein.

        Returns:
        - Protein: The folded protein.

        Raises:
        - ValueError: If a valid folding for the protein cannot be found.
        """
        current = self._protein.get_head()
        self._protein.add_to_grid(current.position, current)

# Start a counter and a movement index to keep track of the movements.
        movement_index: int = 0
        counter: int = 1

        current = current.link

# Loop until a valid folding is found.
        while not self._protein.is_valid():

            # Get the next movement.
            movement = HelixFold.movements[movement_index]

# If the counter is 4, change the movement to a z-movement and
            # reset the counter.
            if counter == 4:
                movement = (0, 0, 1)
                movement_index -= 1
                counter = 0

# Set the position of the next amino acid and add it to the grid.
            current.position = tuple(
                sum(x) for x in zip(current.predecessor.position, movement)
            )
            self._protein.add_to_grid(current.position, current)

# If the next amino acid is the last one, break the loop.
            if current.link is None:
                break

            current = current.link
            movement_index = (movement_index + 1) % len(HelixFold.movements)
            counter += 1

        if not self._protein.is_valid():
            raise ValueError(
                f"Couldn't find a valid folding for protein {self._protein}"
            )

        return self._protein
