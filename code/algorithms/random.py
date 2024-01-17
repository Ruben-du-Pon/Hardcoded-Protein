from typing import Optional, Tuple
from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
import random


class RandomFold:
    """
    Class for performing a random folding on a protein sequence.

    Parameters
    ----------
    protein : Protein
        The protein sequence on which the random folding is performed.
    dimensions : int
        The number of dimensions for the folding (2 for 2D, 3 for 3D).
    no_crossing : Optional[bool], optional
        If True, ensures that no two amino acids occupy the same position,
        by repeatedly selecting a new position until a unique one is found.

    Attributes
    ----------
    _protein : Protein
        The protein sequence on which the random folding is performed.
    _dimensions : int
        The number of dimensions for the folding (2 for 2D, 3 for 3D).
    _no_crossing : Optional[bool]
        If True, ensures that no two amino acids occupy the same position.

    Methods
    -------
    run() -> None:
        Perform the random folding on the protein sequence.
        If no_crossing is True, ensures that no two amino acids occupy the
        same position by repeatedly selecting a new position until a unique one is found.

    get_random_direction() -> Tuple[int, int, int]:
        Get a random direction for the folding.
        Returns a tuple representing the random direction.

    set_position(acid: Aminoacid) -> None:
        Set the position of an amino acid based on a random direction.
        If no_crossing is True, ensures that no two amino acids occupy the same position.

    """  # noqa

    def __init__(self, protein: Protein, dimensions: int,
                 no_crossing: Optional[bool] = False) -> None:
        """
        Initialize the RandomFold object.

        Parameters
        ----------
        protein : Protein
            The protein structure to which the random folding algorithm is applied.
        dimensions : int
            Represents the dimensions of folding (2 or 3).
        no_crossing : Optional[bool], default: False
            If True, ensures that the folding algorithm avoids crossing paths.

        Raises
        ------
        ValueError
            If dimensions is not 2 or 3.

        Notes
        -----
        The RandomFold object is used to apply a random folding algorithm to the given protein.
        The algorithm starts from the second amino acid in the protein sequence and adds each
        amino acid to a grid in a random pattern.

        """
        if dimensions not in (2, 3):
            raise ValueError("Dimensions must be 2 or 3.")

        self._protein = protein
        self._dimensions = dimensions
        self._no_crossing = no_crossing

    def run(self) -> None:
        """
        Perform the random folding on the protein sequence.

        If no_crossing is True, ensures that no two amino acids occupy the
        same position by repeatedly selecting a new position until a unique one is found.
        """  # noqa
        current = self._protein.get_list().link
        while current.link is not None:
            self.set_position(current)
            current = current.link

    def get_random_direction(self) -> Tuple[int, int, int]:
        """
        Get a random direction for the folding.

        Returns
        -------
        Tuple[int, int, int]
            A tuple representing the random direction.
        """
        if self._dimensions == 2:
            directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
        elif self._dimensions == 3:
            directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                          (0, -1, 0), (0, 0, 1), (0, 0, -1)]
        return random.choice(directions)

    def set_position(self, acid: Aminoacid) -> None:
        """
        Set the position of an amino acid based on a random direction.

        If no_crossing is True, ensures that no two amino acids occupy the same position.

        Parameters
        ----------
        acid : Aminoacid
            The amino acid for which the position is set.
        """  # noqa
        while True:
            random_direction = self.get_random_direction()
            new_position = tuple(
                x + y for x, y in zip(acid.position, random_direction))

            if not self._no_crossing or self._protein.is_valid_fold(new_position):
                acid.position = new_position
                self._protein.add_to_grid(new_position, acid)
                break
