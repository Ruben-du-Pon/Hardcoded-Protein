from .aminoacid import Aminoacid
from operator import sub
import csv
from typing import Dict, List, Set, Tuple


class Protein:
    """
    Represents a protein structure composed of amino acids.

    Attributes
    ----------
    _sequence : str
        A string representing the amino acid sequence of the protein.
    _grid : Dict[Tuple[int, int, int], Aminoacid]
        A grid mapping positions to amino acids in the protein structure.
    _head : Aminoacid
        The head of the double-linked list representing the protein structure.
    _score : int
        The stability score of the protein based on its structure.

    Methods
    -------
    __init__(sequence: str) -> None:
        Initialize a Protein object with the given amino acid sequence.

    __create_double_linked_list() -> Aminoacid:
        Create a double-linked list based on the provided amino acid sequence.
        This method is used only during the initialization of the Protein object
        and cannot be called on an already created object.

    get_score() -> int:
        Get the stability score of the protein.

    get_folding() -> List[Dict[str, int]]:
        Get the folding information of the protein.

    create_csv(filename: str) -> None:
        Creates a CSV file that displays a specific folding of a protein.

    is_valid() -> bool:
        Check if every amino acid in the protein has a different position.

    is_valid_fold(position: Tuple[int, int, int]) -> bool:
        Check if a given fold position is valid.

    get_head() -> Aminoacid:
        Get the head of the double-linked list representing the protein structure.

    get_tail() -> Aminoacid:
        Get the tail of the double-linked list representing the protein structure.

    get_grid() -> Dict[Tuple[int, int, int], Aminoacid]:
        Get the protein grid.

    add_to_grid(position: Tuple[int, int, int], acid: Aminoacid) -> None:
        Add an amino acid to the protein grid at the specified position.

    remove_from_grid(position: Tuple[int, int, int]) -> None:
        Remove an amino acid from the protein grid at the specified position.

    reset_grid():
        Clear the grid and add back all the positions of the amino acids.

    __str__() -> str:
        Return the string representation of the protein.

    __len__() -> int:
        Return the length of the amino acid sequence.
    """

    def __init__(self, sequence: str) -> None:
        """
        Initialize a Protein object with the given amino acid sequence.

        Parameters
        ----------
        sequence : str
            A string representing the amino acid sequence of the protein.
        """
        self._sequence: str = sequence
        self._list: List[Aminoacid] = []
        self._grid: Dict[Tuple[int, int, int], Aminoacid] = {}
        self._head: Aminoacid = self.__create_double_linked_list()
        self._score: int = 0

    def __create_double_linked_list(self) -> Aminoacid:
        """
        Create a double-linked list based on the provided amino acid sequence.

        This method is used only during the initialization of the Protein object
        and cannot be called on an already created object.

        Returns
        -------
        Aminoacid
            The head of the double-linked list.
        """
        if not self._sequence:
            return None

        # Create the first amino acid in the sequence
        head = Aminoacid(type=self._sequence[0])
        current = head
        self._list.append(current)

        # Create the rest of the amino acids in the sequence
        for _, amino_type in enumerate(self._sequence[1:], start=1):
            new_aminoacid = Aminoacid(
                type=amino_type, predecessor=current)
            current.link = new_aminoacid
            current = new_aminoacid
            self._list.append(current)

        return head

    def get_score(self) -> int:
        """
        Calculate the stability score of the protein.

        The score is calculated based on the adjacency (non-diagonal)
        of different amino acids in the protein structure. Each adjacent
        pair contributes to the score based on their types.

        Returns
        -------
        int
            The calculated stability score of the protein.
        """
        self._score = 0
        current = self._head

        while current:

            # Get the positions of the amino acids connected to the current amino acid
            connections = [pos for node in [current.predecessor, current.link]
                           if node and (pos := getattr(node, 'position', None))
                           and isinstance(pos, Tuple) and len(pos) == 3]

            # Define the adjacent positions to the current amino acid
            adjacent_positions = [(1, 0, 0), (-1, 0, 0),
                                  (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

            # Check if the adjacent positions are not yet in the connections
            check_positions = [(c + adj[0], d + adj[1], e + adj[2]) for
                               (c, d, e), adj in zip([current.position] *
                                                     len(adjacent_positions),
                                                     adjacent_positions) if
                               (c + adj[0], d + adj[1], e + adj[2]) not in
                               connections]

            # Add the stability score of the current amino acid to the total score
            for pos_tuple in check_positions:
                if pos_tuple in self._grid:
                    self._score += current.get_stability_score(
                        self._grid[pos_tuple])

            current = current.link

        # Return the total score divided by 2 since every connection is counted twice
        return (self._score // 2)

    def get_folding(self) -> List[Dict[str, int]]:
        """
        Get the folding information of the protein.

        Returns
        -------
        List[Dict[str, int]]
            A list of dictionaries containing amino acid types and folding directions,
            with a final entry for the stability score.
        """
        folding = []
        current = self._head

        # Return empty list if the protein is empty
        if current is None or current.link is None:
            return folding

        # Get the folding information of the protein
        while current is not None:

            # Set the fold to 0 if the current amino acid is the last one
            if current.link is None:
                fold = 0

            # Calculate the fold based on the difference between the current
            # and next amino acid
            else:
                next_amino = current.link
                difference = tuple(
                    map(sub, next_amino.position, current.position))
                fold = difference[0] + 2 * difference[1] + 3 * difference[2]

            folding.append({'amino': current.get_type(), 'fold': fold})
            current = current.link

        return folding

    def create_csv(self, filename: str, verbose: bool = False) -> None:
        """
        Creates a CSV file that displays a specific folding of a protein.

        Parameters
        ----------
        filename : str
            The name of the CSV file to be created.

        Postconditions
        --------------
        Creates a CSV file in the "data/output/csv/" directory with folding data.
        The file includes a header "amino", "fold", a footer "score", <score>,
        and a body with amino acid types P, H, or C followed by directions 1, -1, 2, -2, 3, or -3.
        """

        # Create the CSV file
        with open(filename, 'w', newline='') as file:

            # Define the header and writer
            header = ["amino", "fold"]
            writer = csv.DictWriter(file, fieldnames=header)

            # Retrieve the folding information
            folding = self.get_folding()

            # Add the score to the folding information
            folding.append({'amino': 'score', 'fold': self.get_score()})

            # Write the folding information to the CSV file
            writer.writeheader()
            writer.writerows(folding)
            print(f"{filename} created.") if verbose else None

    def is_valid(self) -> bool:
        """
        Check if every amino acid in the protein has a different position and
        the distance between two consecutive aminoacids is always 1.

        Returns
        -------
        bool
            True if every amino acid has a different position and the distance
            between two consecutive aminoacids is always 1, False otherwise.
        """
        current = self._head

        while current.link:

            # Check if the distance between two consecutive aminoacids is always 1
            if abs(current.position[0] - current.link.position[0]) + \
                    abs(current.position[1] - current.link.position[1]) + \
                    abs(current.position[2] - current.link.position[2]) != 1:
                return False
            current = current.link

        # Check if every amino acid has a different position
        return len(self._grid) == len(self._sequence)

    def is_valid_fold(self, position: Tuple[int, int, int]) -> bool:
        """
        Check if a given fold position is valid.

        Parameters
        ----------
        position : Tuple[int, int, int]
            The position to be checked.

        Returns
        -------
        bool
            True if the chosen position is not yet a key in the protein grid, False otherwise.
        """
        return position not in self._grid

    def get_head(self) -> Aminoacid:
        """
        Get the head of the double-linked list representing the protein structure.

        Returns
        -------
        Aminoacid
            The head of the double-linked list.
        """
        return self._head

    def get_tail(self) -> Aminoacid:
        """
        Get the tail of the double-linked list representing the protein structure.

        Returns
        -------
        Aminoacid
            The tail of the double-linked list.
        """
        current = self._head.link
        while current.link is not None:
            current = current.link

        return current

    def get_list(self) -> List[Aminoacid]:
        """
        Get the list of amino acids in the protein.

        Returns
        -------
        List[Aminoacid]
            The list of amino acids in the protein.
        """
        return self._list

    def get_grid(self) -> Dict[Tuple[int, int, int], Aminoacid]:
        """
        Get the protein grid.

        Returns
        -------
        Dict[Tuple[int, int, int], Aminoacid]
            The protein grid mapping positions to amino acids.
        """
        return self._grid

    def add_to_grid(self, position: Tuple[int, int, int], acid: Aminoacid) -> None:
        """
        Add an amino acid to the protein grid at the specified position.

        Parameters
        ----------
        position : Tuple[int, int, int]
            The position where the amino acid should be added.
        acid : Aminoacid
            The amino acid to be added to the grid.
        """
        self._grid[position] = acid

    def remove_from_grid(self, position: Tuple[int, int, int]) -> None:
        """
        Remove an amino acid from the protein grid at the specified position.

        Parameters
        ----------
        position : Tuple[int, int, int]
            The position from which the amino acid should be removed.
        """
        self._grid.pop(position, None)

    def reset_grid(self) -> None:
        """
        Clear the grid and add back all the positions of the amino acids.
        """

        # Clear the grid
        self._grid.clear()

        # Add back all the positions of the amino acids
        current = self.get_head()
        while current:
            self._grid[current.position] = current
            current = current.link

    def __str__(self) -> str:
        """
        Return the string representation of the protein.

        Returns
        -------
        str
            The amino acid sequence of the protein.
        """
        return self._sequence

    def __len__(self) -> int:
        """
        Return the length of the amino acid sequence.

        Returns
        -------
        int
            The length of the amino acid sequence.
        """
        return len(self._sequence)

    def __getstate__(self) -> Tuple[str, List[Aminoacid],
                                    Dict[Tuple[int, int, int], Aminoacid],
                                    Aminoacid, int]:
        """
        Return state values to be pickled.

        Returns
        -------
        Tuple[str, List[Aminoacid], Dict[Tuple[int, int, int], Aminoacid], Aminoacid, int]
            The state values to be pickled.
        """  # noqa
        return (self._sequence, self._list, self._grid,
                self._head, self._score)

    def __setstate__(self, state: Tuple[str, List[Aminoacid],
                                        Dict[Tuple[int, int, int], Aminoacid],
                                        Aminoacid, int]) -> None:
        """
        Restore state from the unpickled state values.

        Parameters
        ----------
        state : Tuple[str, List[Aminoacid], Dict[Tuple[int, int, int], Aminoacid], Aminoacid, int]
            The state values to be unpickled.
        """  # noqa
        self._sequence, self._list, self._grid, self._head, self._score = state
