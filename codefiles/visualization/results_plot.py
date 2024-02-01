import argparse
import os
import matplotlib.pyplot as plt
import csv
import numpy as np


def read_csv_file(file_path):
    with open(file_path, 'r') as file:
        csvreader = csv.reader(file)
        rows = [row for row in csvreader if row]  # Skip empty lines

        number_of_iterations = 1000
        y, avg, names = [], [], []
        seq_values = []

        for row in rows:
            if row[0].startswith('Average:'):
                try:
                    avg_value = float(row[1])
                    avg.append(avg_value)
                    y.append(seq_values)
                    seq_values = []  # Reset for the next sequence
                except ValueError:
                    print(
                        f"Warning: Could not parse average value '{row[1]}' as float.")
                    avg.append(None)
            else:
                try:
                    seq_values.append(float(row[-1]))
                    if len(seq_values) == 1:  # First entry of a new sequence
                        # Assuming the protein chain name is in the second column
                        names.append(row[1])
                except ValueError:
                    print(
                        f"Warning: Could not parse value '{row[-1]}' as float.")

        return y, avg, names


def plot_combined_histogram(data, averages, names, algorithm_names, output_folder, file_name, dimension):
    os.makedirs(output_folder, exist_ok=True)

    # Define a color palette for the algorithms
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd']
    lighter_colors = ['#add8e6', '#ffd699', '#90ee90',
                      '#c5b0d5']  # Lighter shades for 4 algorithms

    for i in range(len(data[0])):
        fig, ax = plt.subplots()

        # Determine the range of stability scores for this iteration
        iteration_scores = [seq_data[i] for seq_data in data]

        min_score, max_score = min(map(min, iteration_scores)), max(
            map(max, iteration_scores))
        max_score = max(0, max_score)  # Ensure the range goes up to at least 0
        score_range = np.arange(min_score, max_score + 1, 1)

        # Width of bars and the number of algorithms
        width = 0.2
        num_algorithms = len(data)

        # Group bars for each stability score
        for j, (seq_data, avg) in enumerate(zip(data, averages)):
            # Count the frequency of each score
            score_counts = [seq_data[i].count(score) for score in score_range]
            ax.bar(score_range - (width * num_algorithms / 2) + j * width,
                   score_counts, width, color=colors[j], label=f"{algorithm_names[j]}")

            # Add a line for the average with a darker color
            ax.axvline(x=avg[i], color=colors[j],
                       linestyle='dashed', linewidth=1.5)

            # Add a line for the minimum with a lighter color
            min_val = min(seq_data[i])
            ax.axvline(
                x=min_val, color=lighter_colors[j], linestyle='dotted', linewidth=1.5)

       # Create custom handles for the legend
        handles, labels = ax.get_legend_handles_labels()
        avg_line = plt.Line2D([0], [0], color='black',
                              linestyle='dashed', linewidth=1.5, label='Avg')
        min_line = plt.Line2D([0], [0], color='gray',
                              linestyle='dotted', linewidth=1.5, label='Min')

        # Insert Avg and Min lines at the beginning of the handles list
        handles = [avg_line, min_line] + handles

        # Set plot title to include the protein sequence name
        protein_name = names[i] if i < len(names) else "Unknown Protein"
        ax.set_title(f"{protein_name} - {dimension}")
        ax.set_xlabel('Stability Score')
        ax.set_ylabel('Frequency')
        ax.legend(handles=handles, fontsize='small',
                  loc='upper left', bbox_to_anchor=(0, 1))
        plt.grid(True)

        fig_name = f"{protein_name}_{file_name}.png"
        fig.savefig(os.path.join(output_folder, fig_name))
        print(f"Saved plot: {os.path.join(output_folder, fig_name)}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate plots for protein chain analysis.")
    parser.add_argument("number", type=int, choices=[
                        2, 3], help="Number representing the dimension (2D or 3D).")
    parser.add_argument("letter", type=str, nargs='?', choices=[
                        'c', 'C', ''], default='', help="Letter representing the type (e.g., C), optional.")

    args = parser.parse_args()

    dimension = args.number
    letter = args.letter.upper()  # Converts 'c' to 'C'

    # File name logic
    if letter:
        file_name = f"{dimension}D_{letter}.csv"
    else:
        file_name = f"{dimension}D.csv"
    folder_paths = ["data/output/fress/", "data/output/bfs/",
                    "data/output/hillclimber_data/", "data/output/annealing_data/"]  # Add your folder paths
    algorithm_names = ["Fress", "BFS", "Hill Climber", "Simulated Annealing"]

    data = []
    averages = []
    for folder in folder_paths:
        file_path = os.path.join(folder, file_name)
        if os.path.exists(file_path):
            y, avg, names = read_csv_file(file_path)
            data.append(y)
            averages.append(avg)
        else:
            print(f"File not found: {file_path}")
            return

    if letter:
        output_path = f'data/output/plot_results/{dimension}D_{letter}'
    else:
        output_path = f'data/output/plot_results/{dimension}D'
    plot_combined_histogram(
        data, averages, names, algorithm_names, output_path, file_name, f"{dimension}D")


if __name__ == "__main__":
    main()
