import matplotlib.pyplot as plt
from typing import Tuple


def plot_3d(protein, colors: Tuple[str, str, str], line_number: int) -> None:
    colors = [color.lower() for color in colors]

    if "C" not in protein._sequence:
        curr_pos = protein._head
        coordinates = [
            (curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])
        ]
        colors_ = [colors[0] if curr_pos._type == "H" else colors[1]]

        while curr_pos.link is not None:
            curr_pos = curr_pos.link
            coordinates.append(
                ((curr_pos.position[0],
                 curr_pos.position[1], curr_pos.position[2]))
            )
            colors_ = colors_ + [colors[0]
                                 if curr_pos._type == "H" else colors[1]]

        x_coordinates = [x for x, _, _ in coordinates]
        y_coordinates = [y for _, y, _ in coordinates]
        z_coordinates = [z for _, _, z in coordinates]

        fig = plt.figure("Protein Alginment")
        ax = fig.add_subplot(projection="3d")

        ax.scatter(x_coordinates, y_coordinates,
                   z_coordinates, s=50, marker="o", c=colors_)

        for i in range(1, protein._length):
            ax.plot(
                [x_coordinates[i - 1], x_coordinates[i]],
                [y_coordinates[i - 1], y_coordinates[i]],
                [z_coordinates[i - 1], z_coordinates[i]],
                color="black",
            )

        # ax.axis("off")
        ax.set_xlim((min(x_coordinates) - 2, max(x_coordinates) + 2))
        ax.set_ylim((min(y_coordinates) - 2, max(y_coordinates) + 2))
        ax.set_zlim((min(z_coordinates) - 2, max(z_coordinates) + 2))
        plt.axis("off")

        legend_labels = ["Hydrofoob", "Polair"]
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color,
                markersize=10,
            )
            for color in colors
        ]
        plt.legend(legend_handles, legend_labels, loc="upper right")

        # plt.show()
        plt.savefig(f"data/output/plot_3D_{line_number}.png")
        print(f"data/output/plot_3D_{line_number}.png created")

    elif "C" in protein._sequence:
        curr_pos = protein._head
        coordinates = [
            (curr_pos.position[0], curr_pos.position[1], curr_pos.position[2])
        ]
        colors_ = [
            colors[0]
            if curr_pos._type == "H"
            else colors[1]
            if curr_pos._type == "P"
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
                if curr_pos._type == "H"
                else colors[1]
                if curr_pos._type == "P"
                else colors[2]
            ]

        x_coordinates = [x for x, _, _ in coordinates]
        y_coordinates = [y for _, y, _ in coordinates]
        z_coordinates = [z for _, _, z in coordinates]

        fig = plt.figure("Protein Alginment")
        ax = fig.add_subplot(projection="3d")

        ax.scatter(x_coordinates, y_coordinates,
                   z_coordinates, s=50, marker="o", c=colors_)

        for i in range(1, protein._length):
            ax.plot(
                [x_coordinates[i - 1], x_coordinates[i]],
                [y_coordinates[i - 1], y_coordinates[i]],
                [z_coordinates[i - 1], z_coordinates[i]],
                color="black",
            )

        ax.set_xlim((min(x_coordinates) - 2, max(x_coordinates) + 2))
        ax.set_ylim((min(y_coordinates) - 2, max(y_coordinates) + 2))
        ax.set_zlim((min(z_coordinates) - 2, max(z_coordinates) + 2))
        plt.axis("off")

        legend_labels = ["Hydrofoob", "Polair", "Cysteine"]
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color,
                markersize=10,
            )
            for color in colors
        ]
        plt.legend(legend_handles, legend_labels, loc="upper right")

        # plt.show()
        plt.savefig(f"data/output/plot_3D_{line_number}.png")
        print(f"data/output/plot_3D_{line_number}.png created")
