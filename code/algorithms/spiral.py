from code.classes.protein import Protein


def spiral_fold(protein: Protein):
    """
    """  # noqa
    current = protein.get_list().link
    protein.add_to_grid(current.position, current)

    movements = [(0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0)]
    movement_index = 0
    # Number of steps in the current direction
    steps = 1

    while True:

        # Flag variable
        loop_break = False

        # Move in the current direction for the specified number of steps
        for _ in range(steps):

            # Set position to the sum of current posi
            current.position = tuple(
                sum(x) for x in zip(current.predecessor.position,
                                    movements[movement_index])
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

    if not protein.is_valid:
        raise ValueError(
            f"Couldn't find a valid folding for protein {protein}")
