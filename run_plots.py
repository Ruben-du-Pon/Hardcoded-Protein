from codefiles.visualization import (
    annealing_plot,
    hillclimber_plot,
    results_plot
)


def main() -> None:
    """
    Main function that runs the program.
    """
    # Run the hillclimber algorithm
    hillclimber_plot.main()

    # Run the annealing algorithm
    annealing_plot.main()

    # Plot the results
    # results_plot.main()


if __name__ == "__main__":
    main()
