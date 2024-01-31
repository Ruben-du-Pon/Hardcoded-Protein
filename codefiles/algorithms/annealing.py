import copy
import random
import csv
from tqdm import tqdm
from typing import List, Optional
from .random import RandomFold
from ..classes.protein import Protein
from .hillclimber import HillclimberFold


class AnnealingFold(HillclimberFold):
    """
    Represents a protein folding algorithm based on Simulated Annealing.

    Inherits from HillclimberFold.

    Attributes:
    - protein (Protein): The protein to fold.
    - dimensions (int): The number of dimensions for the protein fold (2 or 3).
    - iterations (int): The number of iterations for the protein folding process.
    - scores (List[int]): The list of scores obtained during the folding process.
    - outputfile (Optional[str]): The output file to write the folding results.
    - verbose (Optional[bool]): Flag indicating whether to print verbose output.

    Methods:
    - run(): Runs the Simulated Annealing algorithm for protein folding.
    - _check_highscore(protein: Protein) -> bool: Checks if a new protein fold is better than the current highscore.

    Raises:
    - ValueError: If the specified fold algorithm is invalid.
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 scores: List[int] = [],
                 outputfile: Optional[str] = None,
                 verbose: Optional[bool] = False) -> None:
        super().__init__(protein, dimensions, iterations, scores, outputfile,
                         verbose)
        self._temperature = 10.0
        self._cooling_rate = 0.0025

    def run(self) -> Protein:
        """
        Main function that reads protein sequences from a file and generates protein folds using the specified algorithm.

        Usage: python main.py <fold_algorithm> <dimensions> <iterations> [C/c]

        Returns:
        - Protein: The highest scoring protein fold.

        Raises:
        - ValueError: If the specified fold algorithm is invalid.
        """
        # Give feedback that the algorithm has started.
        print("Starting Simulated Annealing fold.")

        # Start with a random fold.
        start_state = RandomFold(self._protein, self._dimensions)
        protein = start_state.run()

        # Start the highscore with the score of the random fold.
        self._highscore = (protein, protein.get_score())

        # Run the algorithm for the specified number of iterations.
        for iteration in tqdm(range(self._iterations)):
            next_fold = copy.deepcopy(self._highscore[0])
            self._run_experiment(next_fold)
            self._temperature *= 1 - self._cooling_rate
            self._temperature = max(1, self._temperature)

            # Write data to file.
            if self._outputfile:

                self._scores.append(next_fold.get_score())

                with open(self._outputfile, "a") as file:
                    writer = csv.writer(file)
                    writer.writerow(
                        [iteration, str(next_fold), next_fold.get_score()])

        # Return the highest scoring protein.
        self._highscore[0].reset_grid()
        return self._highscore[0]

    def _check_highscore(self, protein: Protein) -> bool:
        """
        Checks if the given protein has a higher score than the current highscore protein.

        Parameters:
        - protein (Protein): The protein to be checked.

        Returns:
        - bool: True if the given protein has a higher score and should be accepted as the new highscore, False otherwise.

        Raises:
        - ValueError: If the specified fold algorithm is invalid.
        """
        # Reset the grid.
        protein.reset_grid()

        # Early return if the protein is invalid.
        if not protein.is_valid():
            return False

        # Get the old and new score.
        old_score = self._highscore[1]
        new_score = protein.get_score()

        # Get the temperature and update it.
        temperature = self._temperature

        # Calculate the chance of accepting the new protein.
        chance = 2 ** ((old_score - new_score) / temperature)
        chance = min(2, chance)

        # Check if the new protein is better than the old one and if it should
        # be accepted.
        if (new_score < self._highscore[1]) and (chance > random.random()):
            self._highscore = (protein, new_score)
            return True

        return False
