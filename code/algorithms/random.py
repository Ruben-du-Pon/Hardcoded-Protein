import random
from typing import Optional, Tuple
from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid

DIRECTIONS_2D = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]

DIRECTIONS_3D = [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                 (0, -1, 0), (0, 0, 1), (0, 0, -1)]


class RandomSearch:
    """
    Class for performing a random search on a protein sequence.

    Parameters
    ----------
    protein : Protein
        The protein sequence on which the random search is performed.
    dimensions : int
        The number of dimensions for the search (2 for 2D, 3 for 3D).
    no_crossing : Optional[bool], optional
        If True, ensures that no two amino acids occupy the same position,
        by repeatedly selecting a new position until a unique one is found.

    Attributes
    ----------
    _protein : Protein
        The protein sequence on which the random search is performed.
    _positions : set
        A set containing the positions of amino acids.
    _no_crossing : Optional[bool]
        If True, ensures that no two amino acids occupy the same position.
    _dimensions : int
        The number of dimensions for the search (2 for 2D, 3 for 3D).
    """

    def __init__(self, protein: Protein, dimensions: int,
                 no_crossing: Optional[bool] = False) -> None:
        self._protein = protein
        self._positions: set[Tuple[int, int, int]] = {(0, 0, 0)}
        self._no_crossing: Optional[bool] = no_crossing
        self._dimensions: int = dimensions

    def run(self) -> None:
        """
        Perform the random search on the protein sequence.

        If no_crossing is True, ensures that no two amino acids occupy the
        same position by repeatedly selecting a new position until a unique one is found.
        """
        current = self._protein.get_list().link
        while current.link is not None:
            if self._no_crossing:
                while current.position in self._positions:
                    self.set_position(current)
                self._positions.add(current.position)
                current = current.link
            else:
                self.set_position(current)
                current = current.link

    def get_random_direction(self) -> Tuple[int, int, int]:
        """
        Get a random direction for the search.

        Returns
        -------
        Tuple[int, int, int]
            A tuple representing the random direction.
        """
        if self._dimensions == 2:
            directions = DIRECTIONS_2D
        elif self._dimensions == 3:
            directions = DIRECTIONS_3D
        return random.choice(directions)

    def set_position(self, acid: Aminoacid) -> None:
        """
        Set the position of an amino acid based on a random direction.

        Parameters
        ----------
        acid : Aminoacid
            The amino acid for which the position is set.
        """
        acid.position = tuple(
            x + y + z for x, y, z in zip(acid.predecessor.position,
                                         self.get_random_direction()))
