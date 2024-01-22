from ..classes.protein import Protein
from .random import RandomFold
from typing import List, Tuple, Dict, Optional
import copy
import pandas as pd


class FressFold:
    def __init__(self, protein: Protein, dimensions: int, iterations: int,
                 verbose: Optional[bool] = False) -> None:
        """
        Constructs a FressFold instance for protein folding optimization.
        
        Parameters:
        protein : Protein
            The protein instance on which the folding optimization will be performed.
        """
        self._protein = protein
        self._dimensions = dimensions
        self._verbose = verbose
        self._iterations = iterations

        self.optimization_parameters: Dict[str, float] = {}
        self.stability_scores: List[float] = []

        # keep the best protein fold
        self.best_protein_structure: Protein = protein


    def run(self) -> None:
        """
        Perform the FRESS algortim on a random protein sequence.
        """  # noqa
        best_score = 0
        best_random_score = 0
        best_protein = None
        best_random_protein = None
        for _ in range(5000):
            protein_copy = self._clear_deepcopy()
            fold = RandomFold(protein_copy, self._dimensions, avoid_overlap=True, verbose=False).run()
            current_random_score = fold.get_score()
            if current_random_score < best_random_score:
                best_random_score = current_random_score
                best_random_protein = fold

        print(f"Best random score: {best_random_score}")
        print("Running _analyzeChain()")
        self._analyzeChain(best_random_protein)
        best_protein = best_random_protein # TEMP
        return best_protein

    def _analyzeChain(self, input_protein) -> None:
        """
        Analyzes the protein chain to identify regions for potential optimization.
        This version focuses on the contribution of each individual amino acid.
        """
        print(f"Total score: {input_protein.get_score()}")
        chain_analysis = []
        current = input_protein._head

        while current:
            # Calculate the contribution of the current amino acid
            amino_acid_contribution = 0
            connections = [node for node in [current.predecessor, current.link] if node]

            for conn in connections:
                amino_acid_contribution += current.get_stability_score(conn)

            # Information about the amino acid
            amino_acid_info = {
                'position': current.position,
                'stability_contribution': amino_acid_contribution,
                # Add more details if needed
            }
            chain_analysis.append(amino_acid_info)

            current = current.link

        # Store the analysis in an attribute or process it further
        self.chain_analysis_results = chain_analysis

        # If verbose, print or log the analysis results
        if self._verbose:
            print(pd.DataFrame(chain_analysis)) 

    def _calculate_stability(self) -> float:
        """
        Computes the stability score of the protein's current structure based on its folding.
        
        Returns:
        float
            The calculated stability score for the current protein structure.
        """
        pass

    def _clear_deepcopy(self):
        deepcopy = copy.deepcopy(self._protein)
        current = deepcopy.get_head()
        while current is not None:
            deepcopy.remove_from_grid(current.position)
            current.position = (0, 0, 0)
            current = current.link
        return deepcopy
    
    def optimizeChain(self) -> None:
        """
        Applies optimization techniques to improve the protein's folding and stability.
        """
        pass

    def _get_best_structure(self) -> Protein:
        """
        Retrieves the most stable protein structure discovered during optimization.
        
        Returns:
        Protein
            The protein structure with the best stability score.
        """
        pass

    def _iterate_optimization(self) -> None:
        """
        Conducts a single iteration of the optimization algorithm, adjusting the protein folding as needed.
        """
        pass

    def _update_optimization_parameters(self, new_params: Dict[str, float]) -> None:
        """
        Updates the set of parameters that guide the optimization process.
        
        Parameters:
        new_params : Dict[str, float]
            A dictionary containing the new parameter values to be used in optimization.
        """
        self.optimization_parameters.update(new_params)

