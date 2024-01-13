from code.classes.protein import Protein


def zigzag_fold(protein: Protein):
    """
    """  # noqa
    current = protein.get_list()
    protein.add_to_grid(current.position, current)

    movements = [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, 1, 0)]
    movement_index = 0

    while current.link is not None:
        current = current.link
        current.position = tuple(
            sum(x) for x in zip(current.position,
                                current.predecessor.position,
                                movements[movement_index])
        )
        protein.add_to_grid(current.position, current)

        # Increment movement_index and cycle back to 0 if it exceeds the list length
        movement_index = (movement_index + 1) % len(movements)

    if not protein.is_valid:
        raise ValueError(
            f"Couldn't find a valid folding for protein {protein}")
