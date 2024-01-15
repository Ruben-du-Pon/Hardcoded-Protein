from code.classes.protein import Protein
from typing import List, Tuple


def zigzag_fold(protein: Protein) -> None:
    """
    Apply a zigzag folding algorithm to the given protein.

    The zigzag folding algorithm starts from the first amino acid in the protein sequence,
    and it adds each amino acid to a grid by moving in a zigzag pattern. The algorithm continues
    until the last amino acid in the protein sequence is reached.

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
    zigzag_fold(my_protein)
    """

    current = protein.get_list()
    protein.add_to_grid(current.position, current)

    movements: List[Tuple[int, int, int]] = [
        (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, 1, 0)]
    movement_index: int = 0

    while current.link is not None:
        current = current.link
        current.position = tuple(
            sum(x) for x in zip(current.position, current.predecessor.position, movements[movement_index])
        )
        protein.add_to_grid(current.position, current)

        # Increment movement_index and cycle back to 0 if it exceeds the list length
        movement_index = (movement_index + 1) % len(movements)

    if not protein.is_valid():
        raise ValueError(
            f"Couldn't find a valid folding for protein {protein}")
