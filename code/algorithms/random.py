from typing import List, Optional, Tuple, Set
from operator import add
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
        current = self._protein.get_list()
        self._protein.add_to_grid(current.position, current)
        current = current.link
        while current:
            self.set_position(current, self._no_crossing)
            current = current.link

    def set_position(self, acid: Aminoacid, no_crossing: bool, move_history=None) -> None:
        """
        Set the position of an amino acid based on a random direction.

        If no_crossing is True, ensures that no two amino acids occupy the same position.

        Parameters
        ----------
        acid : Aminoacid
            The amino acid for which the position is set.
        """  # noqa
        if move_history is None:
            move_history = []

        directions = self._get_directions()

        if no_crossing:
            while True:
                valid_directions = [d for d in directions if d not in move_history]
                random.shuffle(valid_directions)

                for direction in valid_directions:
                    if acid.predecessor and acid.predecessor.position:
                        new_position = tuple(map(add, acid.predecessor.position, direction))
                        if self._protein.is_valid_fold(new_position):
                            acid.position = new_position
                            self._protein.add_to_grid(new_position, acid)
                            move_history.append(direction)  # Record successful move
                            return

                if move_history:  # Backtrack if there are moves to undo
                    last_move = move_history.pop()
                    reverse_move = tuple(-x for x in last_move)
                    if acid.position:
                        backtrack_position = tuple(map(add, acid.position, reverse_move))
                        self._protein.remove_from_grid(acid.position)
                        acid.position = backtrack_position
                    if acid.predecessor:
                        self.set_position(acid.predecessor, no_crossing, move_history)
                        return
                else:
                    # No more moves to backtrack or at the start of the chain
                    break
        else:
            if acid.predecessor and acid.predecessor.position:
                direction = random.choice(directions)
                new_position = tuple(map(add, acid.predecessor.position, direction))
                acid.position = new_position
                self._protein.add_to_grid(new_position, acid)


    def _get_directions(self) -> List[Tuple[int, int, int]]:
        if self._dimensions == 2:
            return [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
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
