import matplotlib.pyplot as plt
from typing import Tuple
from ..classes.protein import Protein

def plot_3d(protein: Protein, colors: Tuple[str, str, str], filename: str) -> None:
    """
    Plot a 3D representation of the protein structure.

    The function generates a 3D scatter plot of amino acid positions in the protein structure.
    Amino acids of different types are distinguished by colors. The resulting plot is saved as an image file.

    Parameters
    ----------
    protein : Protein
        The protein structure to be visualized.
    colors : Tuple[str, str, str]
        A tuple of three colors representing different amino acid types (Hydrophobic, Polar, Cysteine).
    filename : str
        The name of the output image file.

    Raises
    ------
    ValueError
        If the provided protein sequence is invalid.

    Notes
    -----
    - The function assumes that the protein structure has a double-linked list representation.
    - Amino acid types are identified by characters "H" (Hydrophobic), "P" (Polar), and "C" (Cysteine).

    Example
    -------
    colors = ("red", "blue", "green")
    protein_sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    my_protein = Protein(protein_sequence)
    plot_3d(my_protein, colors, "my_protein_3d.svg")
    """

    # Clear the previous plot
    plt.clf()

    # Convert colors to lowercase for consistency
    colors = [color.lower() for color in colors]

    if "C" not in protein._sequence:
        # If there are no Cysteine residues in the sequence
        curr_pos = protein.get_head()
        coordinates = [(curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])]
        colors_ = [colors[0] if curr_pos.get_type() == "H" else colors[1]]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(((curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])))
            colors_ = colors_ + [colors[0] if curr_pos.get_type() == "H" else colors[1]]

        x_coordinates = [x for x, _, _ in coordinates]
        y_coordinates = [y for _, y, _ in coordinates]
        z_coordinates = [z for _, _, z in coordinates]

        # Create a new 3D plot
        fig = plt.figure("Protein Alignment")
        ax = fig.add_subplot(projection="3d")

        # Scatter plot for amino acid positions
        ax.scatter(x_coordinates, y_coordinates, z_coordinates, s=50, marker="o", c=colors_)

        # Line plot connecting amino acid positions
        for i in range(1, len(protein)):
            ax.plot([x_coordinates[i - 1], x_coordinates[i]],
                    [y_coordinates[i - 1], y_coordinates[i]],
                    [z_coordinates[i - 1], z_coordinates[i]], color="black")

        # Set plot limits and turn off axis
        ax.set_xlim((min(x_coordinates) - 2, max(x_coordinates) + 2))
        ax.set_ylim((min(y_coordinates) - 2, max(y_coordinates) + 2))
        ax.set_zlim((min(z_coordinates) - 2, max(z_coordinates) + 2))
        plt.axis("off")

        # Create a legend outside plt.legend
        legend_labels = ["H", "P"]  # Replace with your custom characters
        legend_handles = [
            plt.Line2D([0], [0], marker="o", color=color, markerfacecolor=color, markersize=10)
            for color in colors
        ]
        plt.legend(legend_handles, legend_labels, loc="upper right")

        # Add a text annotation for the score
        score_text = f"Score: {protein.get_score()}"
        ax.text(min(x_coordinates) - 2, max(y_coordinates) + 2, max(z_coordinates) + 2,
                score_text, fontsize=12.5, color='red')

        # Save the plot as an SVG file
        plt.savefig(filename, format='svg')
        print(f"{filename} created")

    elif "C" in protein._sequence:
        # If there are Cysteine residues in the sequence
        curr_pos = protein.get_head()
        coordinates = [(curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])]
        colors_ = [colors[0] if curr_pos.get_type() == "H" else colors[1] if curr_pos.get_type() == "P" else colors[2]]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(((curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])))
            colors_ = colors_ + [colors[0] if curr_pos.get_type() == "H" else colors[1] if curr_pos.get_type() == "P" else colors[2]]

        x_coordinates = [x for x, _, _ in coordinates]
        y_coordinates = [y for _, y, _ in coordinates]
        z_coordinates = [z for _, _, z in coordinates]

        # Create a new 3D plot
        fig = plt.figure("Protein Alignment")
        ax = fig.add_subplot(projection="3d")

        # Scatter plot for amino acid positions
        ax.scatter(x_coordinates, y_coordinates, z_coordinates, s=50, marker="o", c=colors_)

        # Line plot connecting amino acid positions
        for i in range(1, len(protein)):
            ax.plot([x_coordinates[i - 1], x_coordinates[i]],
                    [y_coordinates[i - 1], y_coordinates[i]],
                    [z_coordinates[i - 1], z_coordinates[i]], color="black")

        # Set plot limits and turn off axis
        ax.set_xlim((min(x_coordinates) - 2, max(x_coordinates) + 2))
        ax.set_ylim
