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

    calculate_score() -> int:
        Calculate the stability score of the protein based on the adjacency of amino acids.

    get_score() -> int:
        Get the stability score of the protein. Calculates the score if not already calculated.

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

    def calculate_score(self) -> int:
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
        def are_connected(amino1: Aminoacid, amino2: Aminoacid) -> bool:
            """Check if two amino acids are connected in the sequence."""
            return amino1.link == amino2 or amino2.link == amino1

        def calculate_distance(position1: Tuple[int, int, int],
                               position2: Tuple[int, int, int]) -> int:
            """Calculate the Manhattan distance between two positions."""
            return sum(abs(p1 - p2) for p1, p2 in zip(position1, position2))

        visited: Set[Aminoacid] = set()
        current = self._head

        while current:
            if current not in visited:
                visited.add(current)
                for other in visited:
                    if current is not other and not are_connected(current,
                                                                  other):
                        if calculate_distance(current.position,
                                              other.position) == 1:
                            self._score += current.stability_score(other)
            current = current.link

        return self._score

    def get_score(self) -> int:
        """
        Get the stability score of the protein. Calculates the score if not already calculated.

        Returns
        -------
        int
            The stability score of the protein.
        """
        if self._score == 0:
            self.calculate_score()
        return self._score

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

        folding.append({'amino': 'score', 'fold': self.get_score()})
        return folding

    def create_csv(self, index: int = 0, algorithm: str = '') -> None:
        """
        Creates a CSV file that displays a specific folding of a protein.

        Parameters
        ----------
        index : int, optional
            An index used to name the output CSV file. Defaults to 0.
        algorithm : str, optional
            The algorithm used for folding, included in the output filename.

        Postconditions
        --------------
        Creates a CSV file in the "data/output/csv/" directory with folding data.
        The file includes a header "amino", "fold", a footer "score", <score>,
        and a body with amino acid types P, H, or C followed by directions 1, -1, 2, -2, 3, or -3.
        """
        filename = f"data/output/csv/output_{algorithm}_{index}.csv"
        with open(filename, 'w', newline='') as file:
            print(f"{filename} created.")
            header = ["amino", "fold"]
            writer = csv.DictWriter(file, fieldnames=header)
            folding = self.get_folding()

            writer.writeheader()
            writer.writerows(folding)

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

    def get_list(self) -> Aminoacid:
        """
        Get the head of the double-linked list representing the protein structure.

        Returns
        -------
        Aminoacid
            The head of the double-linked list.
        """
        return self._head

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

    def __str__(self) -> str:
        """
        Return the string representation of the protein.

        Returns
        -------
        str
            The amino acid sequence of the protein.
        """
        return self._sequence
