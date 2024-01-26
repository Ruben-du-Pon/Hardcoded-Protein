from ..classes.protein import Protein
from .random import RandomFold
from typing import List, Tuple, Dict, Optional
import copy
import pandas as pd
import math
import random
from operator import add, sub

class FressFold(RandomFold):
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

        self._number_of_improvements = 0


    def run(self) -> None:
        """
        Perform the FRESS algortim on a random protein sequence.
        """  # noqa
        return self._start_random_chain()

    def _start_random_chain(self):
        best_score = 0
        best_random_score = 0
        best_protein = None
        best_random_protein = None
        while best_random_score == 0:
            protein_copy = self._clear_deepcopy()
            fold = RandomFold(protein_copy, self._dimensions, avoid_overlap=True, verbose=False).run()
            current_random_score = fold.get_score()
            if current_random_score < best_random_score:
                best_random_score = current_random_score
                best_random_protein = fold

        print(f"Random start score: {best_random_score}")
        return self._analyzeChain(best_random_protein)


    def _analyzeChain(self, input_protein) -> None:
        """
        Analyzes the protein chain to identify regions for potential optimization.
        This version focuses on the connections of each individual amino acid.
        """
        chain_analysis = []
        current = input_protein._head

        while current:
            connections_info = []
            minus_one_count = 0

            connections = [pos for node in [current.predecessor, current.link]
                        if node and (pos := getattr(node, 'position', None))
                        and isinstance(pos, Tuple) and len(pos) == 3]

            adjacent_positions = [(1, 0, 0), (-1, 0, 0),
                                (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

            check_positions = [(c + adj[0], d + adj[1], e + adj[2]) for
                            (c, d, e), adj in zip([current.position] *
                                                    len(adjacent_positions),
                                                    adjacent_positions) if
                            (c + adj[0], d + adj[1], e + adj[2]) not in
                            connections]

            for pos_tuple in check_positions:
                if pos_tuple in input_protein._grid:
                    stability_score = current.get_stability_score(input_protein._grid[pos_tuple])
                    if stability_score == -1:
                        minus_one_count += 1

                    connection_detail = {
                        'adjacent_position': pos_tuple,
                        'stability_contribution': stability_score
                    }
                    connections_info.append(connection_detail)

            amino_acid_info = {
                'current_position': current.position,
                'connections': connections_info,
                'minus_one_count': minus_one_count
            }
            chain_analysis.append(amino_acid_info)

            current = current.link

        # Convert to DataFrame for easier analysis
        # Note: DataFrame will have nested lists, which may need special handling for some operations
        df_chain_analysis = pd.DataFrame(chain_analysis)

        # Identify 'good points'
        good_points = df_chain_analysis[df_chain_analysis['minus_one_count'] > 0]

        # Store the results
        self.good_points = good_points

        return self._suggest_improvement(df_chain_analysis, input_protein, len(input_protein._sequence))

    def _suggest_improvement(self, df_chain_analysis, input_protein, max_range_length, percentage=0.333):
        max_range_length = len(input_protein._sequence)
        range_length = math.ceil(max_range_length * percentage)

        # Calculate the average number of H's per segment
        total_h_count = input_protein._sequence.count('H')
        avg_h_per_segment = (total_h_count / 3) - 0.1# Since we have 3 segments: START, MIDDLE, END 
        # Define the indices for each segment
        start_index = range(0, range_length)
        middle_index = range(range_length, 2 * range_length)
        end_index = range(max_range_length - range_length, max_range_length)

        # Calculate average stability contribution and count H's for each segment
        start_avg = df_chain_analysis.iloc[start_index]['minus_one_count'].mean()
        middle_avg = df_chain_analysis.iloc[middle_index]['minus_one_count'].mean()
        end_avg = df_chain_analysis.iloc[end_index]['minus_one_count'].mean()

        start_h_count = input_protein._sequence[:range_length].count('H')
        middle_h_count = input_protein._sequence[range_length:2 * range_length].count('H')
        end_h_count = input_protein._sequence[max_range_length - range_length:].count('H')

        # Analyze and suggest improvements
        min_avg = min(start_avg, middle_avg, end_avg)
        #suggestion = "Consider making changes at the MIDDLE of the protein to improve stability."
        if min_avg == start_avg and float(start_h_count) >= avg_h_per_segment:
            suggestion = "Consider making changes at the START of the protein to improve stability."
        elif min_avg == middle_avg and float(middle_h_count) >= avg_h_per_segment:
            suggestion = "Consider making changes in the MIDDLE of the protein to improve stability."
        elif min_avg == end_avg and float(end_h_count) >= avg_h_per_segment:
            suggestion = "Consider making changes at the END of the protein to improve stability."
        else:
            suggestion = " NONE "
        
        if self._verbose:
            print(suggestion)


        if " NONE " not in suggestion and self._number_of_improvements <= 30:
            self._number_of_improvements += 1
            return self._get_refold_range(input_protein, suggestion)
        elif self._number_of_improvements == 0:
            return self._start_random_chain()

        return input_protein

    def _get_refold_range(self, input_protein: Protein, suggestion: str) -> Tuple[int, int]:
        """
        Determines the start and end indices for refolding based on the suggestion.

        Parameters:
        suggestion : str
            The suggestion on where to refold.

        Returns:
        Tuple[int, int]
            The start and end indices for the refolding section.
        """
        print(suggestion)
        if "START" in suggestion:
            refold_range =  (0, math.ceil(len(self._protein._sequence) / 3))
        elif "MIDDLE" in suggestion:
            mid = len(self._protein._sequence) / 2
            refold_range =  (math.floor(mid - len(self._protein._sequence) / 6), math.ceil(mid + len(self._protein._sequence) / 6))
        elif "END" in suggestion:
            refold_range = (math.floor(len(self._protein._sequence) * 2 / 3), len(self._protein._sequence) - 1)

        return self._refold_section(input_protein, refold_range)
        

    def _calculate_stability(self, input_protein: Protein) -> float:
        """
        Computes the stability score of the protein's current structure based on its folding.
        
        Returns:
        float
            The calculated stability score for the current protein structure.
        """
        return input_protein.get_score()
    
    def _refold_section(self, input_protein: Protein, index_range, num_attempts: int = 1000) -> Protein:
        """
        Randomly refolds a section of the protein and keeps the best fold.

        Parameters:
        input_protein : Protein
            The current protein structure.
        num_attempts : int
            The number of refolding attempts to make.

        Returns:
        Protein
            The protein structure with the best fold found in this step.
        """
        best_protein = input_protein
        best_score = self._calculate_stability(input_protein)

        start_index, end_index = index_range  # Determine the refold range based on the suggestion
        #print(start_index, end_index)
        for _ in range(num_attempts):
            protein_copy = copy.deepcopy(best_protein)
            
            # Randomly refold the specified section
            fold = self._random_refold(protein_copy, start_index, end_index)

            # Check if the new fold is better
            new_score = self._calculate_stability(fold)
            if new_score <= best_score:
                print(f"The lowest score found is: {new_score}")
                best_protein = fold
                best_score = new_score

        return self._analyzeChain(best_protein) # REPLACE with self._analyzeChain(best_protein) | best_protein
                                                # For recursion | once

    def _random_refold(self, protein: Protein, start_index: int, end_index: int) -> None:
        # Create a deep copy of the protein to work with
        #protein_copy = self._clear_deepcopy()
        start_acid = self._get_acid_at_index(protein, start_index)
        end_acid = self._get_acid_at_index(protein, end_index)

        # Start and end positions for the refolding range
        start_position = start_acid.predecessor.position if start_acid.predecessor else (0,0,0)
        end_position = end_acid.link.position if end_acid.link else None

        # Perform refold
        #f"{start_index}{start_position}\n{end_index}{end_position}")
        
        index = 0
        acid = protein.get_head()
        ridget_part_front = []
        ridget_part_back = []
        front_directions = []
        back_directions = []
        while acid:
            #print(f"{acid}[{index}]: {acid.position}")
            if index <= start_index:
                ridget_part_front.append(acid.position)
            if index >= end_index:
                ridget_part_back.append(acid.position)
            index += 1
            acid = acid.link


        for i in range(len(ridget_part_front) - 1):
            front_directions.append(tuple(ridget_part_front[i + 1][j] - ridget_part_front[i][j] for j in range(3)))
        for i in range(len(ridget_part_back) - 1):
            back_directions.append(tuple(ridget_part_back[i + 1][j] - ridget_part_back[i][j] for j in range(3)))
        
        middle_directions = self._generate_new_fold(protein, start_index, end_index)

        merged_list = front_directions.copy()
        merged_list.extend(middle_directions)
        merged_list.extend(back_directions)

        all_positions = [(0,0,0)]

        for i in range(0, len(merged_list)):
            all_positions.append(tuple(map(add, all_positions[i], merged_list[i])))
        while len(all_positions) != len(set(all_positions)):
            middle_directions = self._generate_new_fold(protein, start_index, end_index)
            merged_list = front_directions.copy()  
            merged_list.extend(middle_directions)
            merged_list.extend(back_directions)
            all_positions = [(0,0,0)]
            for i in range(0, len(merged_list)):
                all_positions.append(tuple(map(add, all_positions[i], merged_list[i])))
            

        return self._create_new_fold(protein, all_positions)


    def _generate_new_fold(self, protein, start_index, end_index) -> None:
        protein_copy = self._clear_deepcopy()
        acid = protein_copy.get_head()
        index = 0
        backtracking_bool = False
        while index < start_index:
            acid = acid.link
            index +=1 

        path = []
        directions = super()._get_directions()

        merges = 1 if index == 0 or end_index == len(protein_copy._sequence) - 1 else 2
        if index == 0:
            new_position = (0,0,0)
            acid.position = new_position
            path.append(new_position)
            acid = acid.link
            index += 1
        while index < end_index + merges:
            if acid.predecessor:
                while directions:
                    random_direction = random.choice(directions)
                    new_position = tuple(map(add, acid.predecessor.position, random_direction))
                    if self._protein.is_valid_fold(new_position):
                        acid.position = new_position
                        path.append(new_position)
                        break
                    else:
                        # remove direction if not valid
                        if random_direction in directions:
                            directions.remove(random_direction)

                if not directions:
                    # backtracking
                    avoid_direction_path = tuple(map(sub, protein_path[index - 2], protein_path[index - 3]))
                    avoid_direction_failed_path = tuple(map(sub, protein_path[index - 1], protein_path[index - 2]))
                    directions = super()._get_directions()
                    if avoid_direction_path in directions:
                        directions.remove(avoid_direction_path)
                    if avoid_direction_failed_path in directions:
                        directions.remove(avoid_direction_failed_path)
                    print("backtracking needed")
                    index -= 1
                    self._protein.remove_from_grid(acid.position)
                    acid = acid.predecessor
                    backtracking_bool = True

            if not backtracking_bool and acid.predecessor:
                avoid_direction = tuple(-x for x in random_direction)
                directions = super()._get_directions()
                if avoid_direction in directions:
                    directions.remove(avoid_direction)
                index += 1
                acid = acid.link
            backtracking_bool = False
        
        path_directions = []
        for i in range(len(path) - 1):
            path_directions.append(tuple(path[i + 1][j] - path[i][j] for j in range(3)))
        
        return path_directions

    def _create_new_fold(self, protein, path):
        protein_copy = self._clear_deepcopy()
        acid = protein_copy.get_head()

        index = 0 
        while acid:
            new_position = path[index]
            acid.position = new_position
            protein_copy.add_to_grid(new_position, acid)
            #print(f"{acid}[{index}]: {acid.position}") if self._verbose else None
            index += 1
            acid = acid.link

        return protein_copy


    def _get_acid_at_index(self, protein: Protein, index: int):
        acid = protein.get_head()
        current_index = 0
        while acid and current_index < index:
            acid = acid.link
            current_index += 1
        return acid

    def _clear_deepcopy(self):
        deepcopy = copy.deepcopy(self._protein)
        current = deepcopy.get_head()
        while current is not None:
            deepcopy.remove_from_grid(current.position)
            current.position = (0, 0, 0)
            current = current.link
        return deepcopy

