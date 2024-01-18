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
        
        move_history = []
        x = 0
        while current:
            self.set_position(current, self._no_crossing, move_history)
            print(f"{current}[{x}]:{current.position}")

            if move_history and move_history[-1] == 'BACKTRACK':
                move_history.pop()
                current = current.predecessor
                if current is None:
                    break
                x -= 1
            else:
                current = current.link
                x += 1

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
                # Generate a list of valid directions not used yet for the current amino acid
                valid_directions = [d for d in directions if (acid.position, d) not in move_history]
                random.shuffle(valid_directions)

                for direction in valid_directions:
                    if acid.predecessor and acid.predecessor.position:
                        new_position = tuple(map(add, acid.predecessor.position, direction))
                        if self._protein.is_valid_fold(new_position):
                            if acid.position != (0,0,0):  # Remove the old position from the grid
                                self._protein.remove_from_grid(acid.position)
                            acid.position = new_position
                            self._protein.add_to_grid(new_position, acid)
                            move_history.append((acid.position, direction))
                            return

                # Backtracking logic
                if move_history:
                    print("backtrackß")
                    while move_history:
                        last_move = move_history.pop()

                        if last_move == 'BACKTRACK':
                            continue  # Skip backtrack markers

                        backtrack_acid = acid.predecessor
                        x = 0
                        while backtrack_acid is not None and x != 10:
                            # Generate valid directions from the backtrack_acid's position
                            valid_directions = [d for d in directions if (backtrack_acid.position, d) not in move_history]
                            random.shuffle(valid_directions)

                            for direction in valid_directions:
                                backtrack_position = tuple(map(add, backtrack_acid.position, direction))
                                if backtrack_position and self._protein.is_valid_fold(backtrack_position):
                                    # Update the state for successful backtracking
                                    if acid.position != (0,0,0):  # Remove the current position from the grid if set
                                        self._protein.remove_from_grid(acid.position)
                                    acid.position = backtrack_position
                                    self._protein.add_to_grid(backtrack_position, acid)
                                    move_history.append((backtrack_acid.position, direction))  # Record the new move
                                    move_history.append('BACKTRACK')  # Mark thiacid.position and s as a backtrack
                                    return  # Successful backtracking

                            backtrack_acid = backtrack_acid.predecessor  # Move to the next predecessor for further backtracking
                            x += 1
                    # If the loop exits without finding a valid position
                    acid.position = (0, 0, 0)
                    self._protein.add_to_grid((0, 0, 0), acid)
                else:
                    break  # Exit if no moves are left in move_history

        else:
            # Handling the case when no_crossing is False
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
