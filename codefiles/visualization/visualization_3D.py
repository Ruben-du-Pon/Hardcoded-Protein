"""
plot_3d Function
Date of Creation: February 1, 2024
Description: This function generates a 3D scatter plot of amino acid positions in a given protein structure.
             Amino acids of different types (Hydrophobic, Polar, Cysteine) are distinguished by colors.
             The resulting plot is saved as an image file. The function assumes that the protein structure
             has a double-linked list representation.
Developer: Ilyass el Allali
"""

import matplotlib.pyplot as plt
from typing import Tuple
from ..classes.protein import Protein
import numpy as np


def plot_3d(protein: Protein, colors: Tuple[str, str, str], filename: str, output: str="svg") -> None:
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
    plot_3d(my_protein, colors, "my_protein.png")
    """
    # Clear the previous plot
    plt.clf()

    colors = [color.lower() for color in colors]

    if "C" not in protein._sequence:
        curr_pos = protein.get_head()
        coordinates = [
            (curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])
        ]
        colors_ = [colors[0] if curr_pos.get_type() == "H" else colors[1]]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(
                ((curr_pos.position[0],
                 curr_pos.position[1], curr_pos.position[2]))
            )
            colors_ = colors_ + [colors[0]
                                 if curr_pos.get_type() == "H" else colors[1]]

        x_coordinates = [x for x, _, _ in coordinates]
        y_coordinates = [y for _, y, _ in coordinates]
        z_coordinates = [z for _, _, z in coordinates]

        fig = plt.figure("Protein Alignment")
        ax = fig.add_subplot(projection="3d")

        ax.scatter(x_coordinates, y_coordinates,
                   z_coordinates, s=50, marker="o", c=colors_)

        for i in range(1, len(protein)):
            ax.plot(
                [x_coordinates[i - 1], x_coordinates[i]],
                [y_coordinates[i - 1], y_coordinates[i]],
                [z_coordinates[i - 1], z_coordinates[i]],
                color="black",
            )


        check_curr = protein.get_head()
        x_cor, y_cor, z_cor = [], [], []

        for i in range(len(x_coordinates) - 1):

            current = check_curr.link.link
            
            while current is not None:
                diff = tuple(np.array(check_curr.position) - np.array(current.position))
                if check_curr._type == "H" and current._type == "H" and (diff == (1, 0, 0) or diff == (-1, 0, 0) or diff == (0, 1, 0) or diff == (0, -1, 0) or diff == (0, 0, 1) or diff == (0, 0, -1)):
                    coordinates = tuple(np.array(current.position) + np.array(tuple(value / 2 for value in diff)))
                    x_cor.append(coordinates[0])
                    y_cor.append(coordinates[1])
                    z_cor.append(coordinates[2])

                current = current.link
            
            check_curr = check_curr.link

        for x, y, z in zip(x_cor, y_cor, z_cor):
            ax.text(x, y-0.05, z-0.05, '*', fontsize=10, color='black', ha='center', va='center')



        x_min, y_min, z_min = min(x_coordinates), min(
            y_coordinates), min(z_coordinates)
        x_max, y_max, z_max = max(x_coordinates), max(
            y_coordinates), max(z_coordinates)


        ax.set_xlim((x_min - 2, x_max + 2))
        ax.set_ylim((y_min - 2, y_max + 2))
        ax.set_zlim((z_min - 2, z_max + 2))
        plt.axis("off")

        legend_labels = ["H", "P"]  # Replace with your custom characters
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color=color,
                markerfacecolor=color,
                markersize=10,
            )
            for color in colors
        ]
        plt.legend(legend_handles, legend_labels, loc="upper right")

        score_text = f"Score: {protein.get_score()}"

        ax.text(x_min - 2, y_max + 2, z_max + 2,
                score_text, fontsize=12.5, color='red')

        ax.text(0.3, 0.3, 0.4, "(0, 0, 0)", fontsize=4, color='black', ha='right', va='top', fontdict={'fontweight': 'bold', 'style': 'italic'})

        if output == "png":
            plt.savefig(filename, format='png')
        else:
            plt.savefig(filename, format='svg')
        print(f"{filename} created")

    elif "C" in protein._sequence:
        curr_pos = protein.get_head()
        coordinates = [
            (curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])
        ]
        colors_ = [
            colors[0]
            if curr_pos.get_type() == "H"
            else colors[1]
            if curr_pos.get_type() == "P"
            else colors[2]
        ]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(
                ((curr_pos.position[0],
                 curr_pos.position[1], curr_pos.position[2]))
            )
            colors_ = colors_ + [
                colors[0]
                if curr_pos.get_type() == "H"
                else colors[1]
                if curr_pos.get_type() == "P"
                else colors[2]
            ]

        x_coordinates = [x for x, _, _ in coordinates]
        y_coordinates = [y for _, y, _ in coordinates]
        z_coordinates = [z for _, _, z in coordinates]

        fig = plt.figure("Protein Alignment")
        ax = fig.add_subplot(projection="3d")

        ax.scatter(x_coordinates, y_coordinates,
                   z_coordinates, s=50, marker="o", c=colors_)

        for i in range(1, len(protein)):
            ax.plot(
                [x_coordinates[i - 1], x_coordinates[i]],
                [y_coordinates[i - 1], y_coordinates[i]],
                [z_coordinates[i - 1], z_coordinates[i]],
                color="black",
            )

        check_curr = protein.get_head()
        x_cor_H, x_cor_C, y_cor_H, y_cor_C, z_cor_H, z_cor_C = [], [], [], [], [], []

        for i in range(len(x_coordinates) - 1):

            current = check_curr.link.link
            
            while current is not None:
                diff = tuple(np.array(check_curr.position) - np.array(current.position))
                if (((check_curr._type == "H" and current._type == "H") or (check_curr._type == "C" and current._type == "H") or (check_curr._type == "H" and current._type == "C")) and (diff == (1, 0, 0) or diff == (-1, 0, 0) or diff == (0, 1, 0) or diff == (0, -1, 0))):
                    coordinates = tuple(np.array(current.position) + np.array(tuple(value / 2 for value in diff)))
                    x_cor_H.append(coordinates[0])
                    y_cor_H.append(coordinates[1])
                    z_cor_H.append(coordinates[2])
                elif check_curr._type == "C" and current._type == "C" and (diff == (1, 0, 0) or diff == (-1, 0, 0) or diff == (0, 1, 0) or diff == (0, -1, 0)):
                    coordinates = tuple(np.array(current.position) + np.array(tuple(value / 2 for value in diff)))
                    x_cor_C.append(coordinates[0])
                    y_cor_C.append(coordinates[1])
                    z_cor_C.append(coordinates[2])

                current = current.link
            
            check_curr = check_curr.link

        for x, y, z in zip(x_cor_H, y_cor_H, z_cor_H):
            ax.text(x, y-0.05, z-0.05, '*', fontsize=12, color='black', ha='center', va='center')


        for z, y, z in zip(x_cor_C, y_cor_C, z_cor_C):
            ax.text(x, y-0.05, z-0.05, '#', fontsize=9, color='black', ha='center', va='center')

        x_min, y_min, z_min = min(x_coordinates), min(
            y_coordinates), min(z_coordinates)
        x_max, y_max, z_max = max(x_coordinates), max(
            y_coordinates), max(z_coordinates)


        ax.set_xlim((x_min - 2, x_max + 2))
        ax.set_ylim((y_min - 2, y_max + 2))
        ax.set_zlim((z_min - 2, z_max + 2))
        plt.axis("off")

        legend_labels = ["H", "P", "C"]  # Replace with your custom characters
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color=color,
                markerfacecolor=color,
                markersize=10,
            )
            for color in colors
        ]
        plt.legend(legend_handles, legend_labels, loc="upper right")

        score_text = f"Score: {protein.get_score()}"
        ax.text(x_min - 2, y_max + 2, z_max + 2,
                score_text, fontsize=12.5, color='red')

        if output == "png":
            plt.savefig(filename, format='png')
        else:
            plt.savefig(filename, format='svg')
        print(f"{filename} created")
