from typing import List, Optional, Tuple, Set
from ..classes.protein import Protein
from ..classes.aminoacid import Aminoacid
import random
import numpy as np


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
                 no_crossing: Optional[bool] = True) -> None:
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

        return self.set_position(self._protein, self._protein._head)



    def set_position(self, prt: Protein, acid: Aminoacid) -> None:
        directions = self._get_directions()
        acid.position = (0, 0, 0)
        current = acid.link
        while current is not None:

            current.position = tuple(np.array(random.choice(directions)) + np.array(current.predecessor.position))
            directions = self._get_directions()

            directions = self.remove_possible_directions(current, current.predecessor, directions)
            current = current.link

        if self.check_valid_fold(prt) is False:
            current = prt.get_list()

            while current is not None:
                prt.remove_from_grid(current.position)
                current = current.link

            return self.set_position(prt, prt._head)
        
        return prt


    def remove_possible_directions(self, aminoacid_curr, aminoacid_pred, directions):

        if sum(np.array(aminoacid_curr.position)) < sum(np.array(aminoacid_pred.position)):
            cant_take_direction = tuple(-(np.array(aminoacid_curr.position) - np.array(aminoacid_pred.position)))
        else:
            cant_take_direction = tuple((np.array(aminoacid_pred.position) - np.array(aminoacid_curr.position)))

        index = directions.index(cant_take_direction)
        del directions[index]
        return directions


    def check_valid_fold(self, prt: Protein):
        current = prt.get_list()
        while current is not None:
            prt.add_to_grid(current.position, current)
            current = current.link

        return prt.is_valid()


    # def backtracking(self, prt: Protein, tail: Aminoacid, directions):
    #     if len(directions) == 0:
    #         directions = self._get_directions()
    #         tail = tail.predecessor
    #         prt.remove_from_grid(tail.position)
    #         tail.position


    #     prt.remove_from_grid(tail.position)

    #     tail.position = tuple(np.array(random.choice(directions)) + np.array(tail.predecessor.position))
    #     directions = self.remove_possible_directions(self, tail, tail.predecessor, directions)
    #     if prt.is_valid():
    #         return prt
    #     else:
    #         return self.backtracking(prt, tail, directions)


    def _get_directions(self) -> List[Tuple[int, int, int]]:
        if self._dimensions == 2:
            return [(1, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0)]
        elif self._dimensions == 3:
            return [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                    (0, -1, 0), (0, 0, 1), (0, 0, -1)]
        return []


    def get_random_direction(self, directions: List[Tuple[int, int, int]]) -> Tuple[int, int, int]:
        """
        Get a random direction for the folding.

        Returns
        -------
        Tuple[int, int, int]
            A tuple representing the random direction.
        """
        return random.choice(directions)
