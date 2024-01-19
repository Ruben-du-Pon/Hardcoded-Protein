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

    create_csv(index: int = 0, algorithm: str = '') -> None:
        Creates a CSV file that displays a specific folding of a protein.

    is_valid() -> bool:
        Check if every amino acid in the protein has a different position.

    is_valid_fold(position: Tuple[int, int, int]) -> bool:
        Check if a given fold position is valid.

    get_list() -> Aminoacid:
        Get the head of the double-linked list representing the protein structure.

    add_to_grid(position: Tuple[int, int, int], acid: Aminoacid) -> None:
        Add an amino acid to the protein grid at the specified position.

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
        self._grid: Dict[Tuple[int, int, int], Aminoacid] = {}
        self._head: Aminoacid = self.__create_double_linked_list()
        self._tail: Aminoacid = self.get_tail()
        self._score: int = 0

    def __len__(self):
        return len(self._sequence)

    def __create_double_linked_list(self):
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

        head = Aminoacid(type=self._sequence[0])
        current = head

        for _, type in enumerate(self._sequence[1:], start=1):
            new_aminoacid = Aminoacid(
                type=type, predecessor=current)
            current.link = new_aminoacid
            current = new_aminoacid

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
            # Add current.predecessor.position and current.link.position to the
            # connections list if they are not None
            connections = [pos for node in [current.predecessor, current.link]
                           if node and (pos := getattr(node, 'position', None))
                           and isinstance(pos, Tuple) and len(pos) == 3]

            adjacent_positions = [(1, 0, 0), (-1, 0, 0),
                                  (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

            # Calculate all positions adjacent to current.position and subtract
            # the positions in connections
            check_positions = [(c + adj[0], d + adj[1], e + adj[2]) for
                               (c, d, e), adj in zip([current.position] *
                                                     len(adjacent_positions),
                                                     adjacent_positions) if
                               (c + adj[0], d + adj[1], e + adj[2]) not in
                               connections]

            for tuple in check_positions:
                if tuple in self._grid:
                    self._score += current.stability_score(self._grid[tuple])

            current = current.link

        return (self._score // 2)

        # def are_connected(amino1: Aminoacid, amino2: Aminoacid) -> bool:
        #     """Check if two amino acids are connected in the sequence."""
        #     return amino1.link == amino2 or amino2.link == amino1

        # def calculate_distance(position1: Tuple[int, int, int],
        #                        position2: Tuple[int, int, int]) -> int:
        #     """Calculate the Manhattan distance between two positions."""
        #     return sum(abs(p1 - p2) for p1, p2 in zip(position1, position2))

        # visited: Set[Aminoacid] = set()
        # current = self._head

        # while current:
        #     if current not in visited:
        #         visited.add(current)
        #         for other in visited:
        #             if current is not other and not are_connected(current,
        #                                                           other):
        #                 if calculate_distance(current.position,
        #                                       other.position) == 1:
        #                     self._score += current.stability_score(other)
        #     current = current.link

        # return self._score

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

        if current is None or current.link is None:
            return folding  # Return empty list in this case

        while current is not None:
            # If current is the last amino acid, there's no next position
            if current.link is None:
                fold = 0
            else:
                next_amino = current.link
                difference = tuple(
                    map(sub, next_amino.position, current.position))
                fold = difference[0] + 2 * difference[1] + 3 * difference[2]

            folding.append({'amino': current.get_type(), 'fold': fold})
            current = current.link
        return folding

    def create_csv(self, filename: str) -> None:
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
        with open(filename, 'w', newline='') as file:
            header = ["amino", "fold"]
            writer = csv.DictWriter(file, fieldnames=header)
            folding = self.get_folding()
            folding.append({'amino': 'score', 'fold': self.get_score()})

            writer.writeheader()
            writer.writerows(folding)
            print(f"{filename} created.")

    def is_valid(self) -> bool:
        """
        Check if every amino acid in the protein has a different position.

        Returns
        -------
        bool
            True if every amino acid has a different position, False otherwise.
        """
        return len(self._sequence) == len(self._grid)

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
        current = self._head.link
        while current.link is not None:
            current = current.link

        return current

    def get_grid(self) -> Dict[Tuple[int, int, int], Aminoacid]:
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
        if position in self._grid:
            # Remove the amino acid from the grid and update its position to None
            removed_acid = self._grid.pop(position)
            removed_acid.position = None

    def __str__(self) -> str:
        """
        Return the string representation of the protein.

        Returns
        -------
        str
            The amino acid sequence of the protein.
        """
        return self._sequence
