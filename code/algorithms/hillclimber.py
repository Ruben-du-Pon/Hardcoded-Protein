import copy
import random
from typing import List, Optional, Tuple
from .random import RandomFold
from ..classes.aminoacid import Aminoacid
from ..classes.protein import Protein
from ..visualization import visualization_2D  # , visualization_3D


class HillclimberFold:
    """
    """

    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 verbose: Optional[bool] = False) -> None:
        """
        """
        if dimensions not in (2, 3):
            raise ValueError(
                "Please enter dimensions as 2 or 3 for 2D or 3D folding.")

        self._protein = protein
        self._dimensions = dimensions
        self._iterations = iterations
        self._highscore = (None, None)
        self._verbose = verbose

    def run(self) -> Protein:
        """
        """
        # Start with a random fold
        start_state = RandomFold(self._protein, self._dimensions, True)
        protein = start_state.run()

        visualization_2D.plot_2d(
            protein, ("red", "blue", "green"), "data/output/plot/hillclimber_start.png", "png")

        protein.create_csv("data/output/csv/hillclimber_start.csv")

        # Start the highscore with the score of the random fold
        starting_fold = copy.deepcopy(protein)
        self._highscore = (starting_fold, starting_fold.get_score())

        # Run the algorithm for the specified number of iterations
        for _ in range(self._iterations):
            print(f"Iteration {_ + 1}") if self._verbose else None
            protein = self._run_experiment(protein)

        visualization_2D.plot_2d(
            protein, ("red", "blue", "green"), "data/output/plot/hillclimber_end.png", "png")

        protein.create_csv("data/output/csv/hillclimber_end.csv")

        print(f"Final highscore: {self._highscore}")
        visualization_2D.plot_2d(
            self._highscore[0], ("red", "blue", "green"), "data/output/plot/hillclimber_highscore.png", "png")

        print(self._highscore[0].is_valid()) if self._verbose else None
        positions = []
        current = self._highscore[0].get_head()
        while current:
            positions.append(current.position)
            current = current.link
        print(positions) if self._verbose else None
        # Return the highest scoring protein
        return self._highscore[0]

    def _run_experiment(self, protein: Protein) -> Protein:
        """
        """
        # Get a random snippet of the protein and give it an intial fold
        snippet, start_position = self._get_snippet(protein)
        protein_copy = copy.deepcopy(protein)
        snippet, changed_protein = self._fold_snippet(
            snippet, start_position, protein_copy)

        # Fold the snippet until a valid fold is found or 100 attemps were made
        count = 0
        while count < 50 and not changed_protein.is_valid():
            snippet, changed_protein = self._fold_snippet(
                snippet, start_position, changed_protein)
            count += 1

        # Change the protein to the new fold if it is a new highscore
        if changed_protein.is_valid() and self._check_highscore(changed_protein):
            return self._highscore[0]

        return protein

    def _get_snippet(self, protein: Protein) -> Tuple[List[Aminoacid], int]:
        """
        """
        # Initialize snippet
        snippet = []

        # Get random length and start position
        length = random.randint(3, len(protein))
        length = min(length, 20)
        start_position = random.randint(0, len(protein) - length)

        current = protein.get_head()

        # Find the start position
        for _ in range(start_position):
            current = current.link

        # Add the amino acids to the snippet
        for _ in range(length):
            snippet.append(current)
            current = current.link

        while current:
            current = current.link

        return snippet, start_position

    def _fold_snippet(self, snippet: List[Aminoacid], snippet_start: int, protein: Protein) -> Tuple[List[Aminoacid], Protein]:
        """
        """
        # Copy the snippet and protein
        snippet_copy = copy.deepcopy(snippet)
        protein_copy = copy.deepcopy(protein)

        # Remove the first amino acid from the snippet
        snippet.pop(0)

        # Create a protein object from the snippet
        sequence = ''.join(str(acid) for acid in snippet)
        snippet_protein = Protein(sequence)

        # Fold the snippet randomly
        RandomFold(snippet_protein, self._dimensions, True).run()

        # Define the translations
        translations = self._find_translations(
            snippet_copy[0], snippet_protein.get_head())

        # Translate the positions of the acids in snippet_protein
        for translation in translations:
            current = snippet_protein.get_head()
            while current:
                current.position = tuple(
                    map(sum, zip(current.position, translation)))
                current = current.link

            # Check if the protein is still a valid fold
            if snippet_protein.is_valid():
                # Set the snippet positions to the new positions
                current = snippet_protein.get_head()
                x = 0
                while current:
                    snippet[x].position = current.position
                    x += 1
                    current = current.link

                # Translate the positions of the rest of the acids in the protein
                new_current = protein.get_head()
                for _ in range(snippet_start):
                    new_current = new_current.link

                for index in range(len(snippet)):
                    new_current.position = snippet[index].position
                    new_current = new_current.link

                second_translations = self._find_translations(
                    new_current, snippet_protein.get_tail())
                for second_translation in second_translations:
                    while new_current:
                        new_current.position = tuple(
                            map(sum, zip(new_current.position, second_translation)))
                        new_current = new_current.link

                    if protein.is_valid():
                        print("Valid folding found!")
                        return snippet, protein

        return snippet_copy, protein_copy

    def _find_translations(self, acid1: Aminoacid,
                           acid2: Aminoacid) -> List[Tuple[int, int, int]]:
        """
        """
        directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0) if
                      self._dimensions == 2 else (1, 0, 0), (-1, 0, 0),
                      (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

        positions = [tuple(x + y for x, y in zip(acid1.position, direction))
                     for direction in directions]

        translations = [tuple(x - y for x, y in zip(acid2.position, position))
                        for position in positions]

        return translations

    def _check_highscore(self, protein: Protein) -> bool:
        """
        """
        # Reset the grid
        protein.reset_grid()

        # Print the current highscore and score if verbose is True
        if self._verbose:
            print(f"Highscore: {self._highscore}")
            print(f"Current highscore: {self._highscore[1]}")
            print(f"Current score: {protein.get_score()}")

        # Check if the protein is a new highscore
        if protein.get_score() < self._highscore[1]:
            self._highscore = (protein, protein.get_score())
            if self._verbose:
                print(f"New highscore found: {self._highscore[1]}")
            return True

        return False
