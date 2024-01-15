from code.classes.protein import Protein
from typing import List, Tuple


def spiral_fold(protein: Protein) -> None:
    """
    Apply a spiral folding algorithm to the given protein.

    The spiral folding algorithm starts from the second amino acid in the protein sequence,
    and it adds each amino acid to a grid in a spiral pattern. The algorithm continues until
    every amino acid in the protein has a different position on the grid.

    Parameters
    ----------
    protein : Protein
        The protein structure to which the folding algorithm is applied.

    Raises
    ------
    ValueError
        If a valid folding cannot be found for the given protein.

    Notes
    -----
    The function assumes that the protein structure has a double-linked list representation,
    and the positions of amino acids are stored as tuples of three integers.

    Example
    -------
    protein_sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    my_protein = Protein(protein_sequence)
    spiral_fold(my_protein)
    """

    current = protein.get_list().link
    protein.add_to_grid(current.position, current)

    movements: List[Tuple[int, int, int]] = [
        (0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0)]
    movement_index: int = 0
    # Number of steps in the current direction
    steps: int = 1

    while not protein.is_valid():

        # Flag variable
        loop_break: bool = False

        # Move in the current direction for the specified number of steps
        for _ in range(steps):

            # Set position to the sum of the current position and the movement vector
            current.position = tuple(
                sum(x) for x in zip(current.predecessor.position, movements[movement_index])
            )
            protein.add_to_grid(current.position, current)
            if current.link is None:
                loop_break = True
                break

            current = current.link

        if loop_break:
            break

        # Increment movement_index and cycle back to 0 if it exceeds the list length
        movement_index = (movement_index + 1) % len(movements)

        # Increase the number of steps every two directions
        if movement_index % 2 == 0:
            steps += 1

    if not protein.is_valid():
        raise ValueError(
            f"Couldn't find a valid folding for protein {protein}")
