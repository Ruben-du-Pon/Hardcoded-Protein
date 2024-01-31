import copy
import csv
import random
from typing import List, Optional, Tuple
# from .random import RandomFold
from .spiral import SpiralFold
from .bfs import BfsFold
from ..classes.protein import Protein
from tqdm import tqdm


class HillclimberFold:
    """
    Represents a hillclimber folding algorithm for proteins.

    Attributes:
    - protein (Protein): The protein to fold.
    - dimensions (int): The number of dimensions for folding (2 or 3).
    - iterations (int): The number of iterations to run the algorithm.
    - scores (List[int]): The list of scores obtained during the algorithm.
    - outputfile (Optional[str]): The path to the output file to write the scores.
    - verbose (Optional[bool]): Whether to print verbose output during the algorithm.

    Methods:
    - run(): Runs the hillclimber folding algorithm and returns the highest scoring protein fold.
    - get_scores(): Returns the list of scores obtained during the algorithm.
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 scores: List[int] = [],
                 outputfile: Optional[str] = None,
                 verbose: Optional[bool] = False) -> None:
        """
        Initializes a new instance of the HillclimberFold class.

        Parameters:
        - protein (Protein): The protein to fold.
        - dimensions (int): The number of dimensions for folding (2 or 3).
        - iterations (int): The number of iterations to run the algorithm.
        - scores (List[int]): The list of scores obtained during the algorithm.
        - outputfile (Optional[str]): The path to the output file to write the scores.
        - verbose (Optional[bool]): Whether to print verbose output during the algorithm.

        Raises:
        - ValueError: If the dimensions parameter is not 2 or 3.
        """
        if dimensions not in (2, 3):
            raise ValueError(
                "Please enter dimensions as 2 or 3 for 2D or 3D folding.")

        self._protein = protein
        self._dimensions = dimensions
        self._iterations = iterations
        self._highscore = (protein, 0)
        self._scores = scores
        self._outputfile = outputfile
        self._verbose = verbose

    def run(self) -> Protein:
        """
        Runs the hillclimber folding algorithm and returns the highest scoring protein fold.

        Returns:
        - Protein: The highest scoring protein fold.

        Raises:
        - ValueError: If the specified fold algorithm is invalid.
        """
        # Implementation details omitted for brevity
        pass

    def get_scores(self) -> List[int]:
        """
        Returns the list of scores obtained during the algorithm.

        Returns:
        - List[int]: The list of scores obtained during the algorithm.
        """
        return self._scores


class HillclimberFold:
    """
    Represents a HillClimber folding algorithm for proteins.

    Attributes:
    - protein (Protein): The protein object to fold.
    - dimensions (int): The number of dimensions for folding (2 or 3).
    - iterations (int): The number of iterations to perform.
    - scores (List[int]): List of scores (optional).
    - outputfile (Optional[str]): Output file path (optional).
    - verbose (Optional[bool]): Whether to print verbose output (optional).

    Methods:
    - __init__(self, protein: Protein, dimensions: int, iterations: int,
               scores: List[int] = [], outputfile: Optional[str] = None,
               verbose: Optional[bool] = False) -> None:
        Initializes a HillClimberFold object.
    - run(self) -> Protein:
        Main function that runs the HillClimber folding algorithm.
    - _run_experiment(self, protein: Protein) -> None:
        Runs an experiment for a given protein.
    - _get_snippet(self, protein: Protein) -> Tuple[int, int]:
        Gets a random snippet of the protein.
    - _process_snippet(self, args: Tuple[Protein, Protein,
                                         Tuple[int, int, int],
                                         Tuple[int, int, int], int]) -> Protein:
        Processes a snippet of the protein.
    - _check_highscore(self, protein: Protein) -> bool:
        Checks if the given protein is a new highscore.
    - get_scores(self) -> List[int]:
        Returns the list of scores.

    Raises:
    - ValueError: If dimensions is not 2 or 3.
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 scores: List[int] = [], outputfile: Optional[str] = None,
                 verbose: Optional[bool] = False) -> None:
        """
        Initialize a HillClimberFold object.

        Parameters:
        - protein (Protein): The protein object to fold.
        - dimensions (int): The number of dimensions for folding (2 or 3).
        - iterations (int): The number of iterations to perform.
        - scores (List[int]): List of scores (optional).
        - outputfile (Optional[str]): Output file path (optional).
        - verbose (Optional[bool]): Whether to print verbose output (optional).

        Raises:
        - ValueError: If dimensions is not 2 or 3.

        Returns:
        - None
        """
        if dimensions not in (2, 3):
            raise ValueError(
                "Please enter dimensions as 2 or 3 for 2D or 3D folding.")

        self._protein = protein
        self._dimensions = dimensions
        self._iterations = iterations
        self._highscore = (protein, 0)
        self._scores = scores
        self._outputfile = outputfile
        self._verbose = verbose

    def run(self) -> Protein:
        """
        Main function that runs the HillClimber folding algorithm.

        Returns:
        - Protein: The highest scoring protein fold.

        Raises:
        - ValueError: If the specified fold algorithm is invalid.
        """
        # Implementation details omitted for brevity
        pass

    def _run_experiment(self, protein: Protein) -> None:
        """
        Runs an experiment for a given protein.

        Parameters:
        - protein (Protein): The protein to run the experiment on.

        Returns:
        - None
        """
        # Implementation details omitted for brevity
        pass

    def _get_snippet(self, protein: Protein) -> Tuple[int, int]:
        """
        Gets a random snippet of the protein.

        Parameters:
        - protein (Protein): The protein to get the snippet from.

        Returns:
        - Tuple[int, int]: The start and end positions of the snippet.
        """
        # Implementation details omitted for brevity
        pass

    def _process_snippet(self, args: Tuple[Protein, Protein,
                                           Tuple[int, int, int],
                                           Tuple[int, int, int], int]) -> Protein:
        """
        Processes a snippet of the protein.

        Parameters:
        - args (Tuple[Protein, Protein, Tuple[int, int, int], Tuple[int, int, int], int]):
          The arguments required for processing the snippet.

        Returns:
        - Protein: The best protein fold obtained from processing the snippet.
        """
        # Implementation details omitted for brevity
        pass

    def _check_highscore(self, protein: Protein) -> bool:
        """
        Checks if the given protein is a new highscore.

        Parameters:
        - protein (Protein): The protein to check.

        Returns:
        - bool: True if the protein is a new highscore, False otherwise.
        """
        # Implementation details omitted for brevity
        pass

    def get_scores(self) -> List[int]:
        """
        Returns the list of scores.

        Returns:
        - List[int]: The list of scores.
        """
        return self._scores
