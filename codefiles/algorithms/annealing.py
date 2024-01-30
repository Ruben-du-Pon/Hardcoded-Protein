import copy
import random
import csv
from tqdm import tqdm
from typing import List, Optional
# from .random import RandomFold
from .spiral import SpiralFold
from ..classes.protein import Protein
from .hillclimber import HillclimberFold


class AnnealingFold(HillclimberFold):
    """
    A class that folds a protein using a simulated annealing algorithm.

    Attributes
    ----------
    _protein : Protein
        The protein to fold.
    _dimensions : int
        The dimensions of the protein.
    _iterations : int
        The number of iterations to run the algorithm.
    _highscore : Tuple[Protein, int]
        The highest scoring protein and its score.
    _verbose : bool
        If True, print additional output.
    _temperature : float
        The initial temperature for the annealing process.
    _cooling_rate : float
        The rate at which the temperature decreases.

    Methods
    -------
    run()
        Runs the simulated annealing algorithm.
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 scores: List[int] = [],
                 outputfile: Optional[str] = None,
                 verbose: Optional[bool] = False) -> None:
        """
        Initializes an AnnealingFold object.

        Parameters
        ----------
        protein : Protein
            The protein to fold.
        dimensions : int
            The dimensions of the protein.
        iterations : int
            The number of iterations to run the algorithm.
        verbose : bool, optional
            If True, print additional output.

        Raises
        ------
        ValueError
            If the dimensions are not 2 or 3.
        """
        super().__init__(protein, dimensions, iterations, scores, outputfile,
                         verbose)
        self._temperature = 10.0
        self._cooling_rate = 0.0025

    def run(self) -> Protein:
        """
        Runs the Simulated Annealing fold algorithm.

        Returns
        -------
        Protein
            The highest scoring protein found.
        """
        # Give feedback that the algorithm has started.
        print("Starting Simulated Annealing fold.")

        # Start with a random fold.
        start_state = SpiralFold(self._protein, self._dimensions)
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
        Checks if the current protein's score is a new highscore.

        Parameters
        ----------
        protein : Protein
            The protein to check.

        Returns
        -------
        bool
            True if the current protein's score is a new highscore, False otherwise.
        """  # noqa
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
        if (new_score < self._highscore[1]) \
                and (chance > random.random()):
            self._highscore = (protein, new_score)
            return True

        return False
