import matplotlib.pyplot as plt
import csv
import math
import numpy as np


def reading_the_csv_file(file):
    """
    Read data from a CSV file and organize it into lists.

    Parameters:
    - file (File): A file object containing CSV data.

    Returns:
    - c (list): List of unique values from the second-to-last column of each sequence.
    - y (list): List of values from the last column of each sequence.
    - avg (list): List of average values from each sequence.
    - amount_seq (int): Number of sequences in the file.
    """
    csvreader = csv.reader(file)

    rows = []
    for row in csvreader:
        rows.append(row)

    number_of_gens = 10**5
    x = number_of_gens
    c, y, avg, amount_seq = [], [], [], math.floor(
        len(rows) / (number_of_gens - 1))

    for _ in range(amount_seq):
        for i in rows[x - number_of_gens:x]:
            y.append(int(i[-1]))
            if i[-2] not in c:
                c.append(i[-2])
        avg.append(rows[x + 1])
        x += number_of_gens + 4

    return c, y, avg, amount_seq


def plot_random(file):
    """
    Plot histograms for random sequences from a CSV file.

    Parameters:
    - file (File): A file object containing CSV data.

    Returns:
    - lst (list): List of matplotlib figure objects.
    """
    c, y, avg, amount_seq = reading_the_csv_file(file)
    lst = []

    number_of_gens = 10**5

    for count, val in enumerate(range(number_of_gens, (amount_seq + 1) * number_of_gens, number_of_gens)):
        # Create a new Figure and Axes for each plot
        fig, ax = plt.subplots()

        # Start bins at the lesser of the minimum value or 0
        bin_start = min(min(y[val - number_of_gens:val]), 0)
        # End bins just beyond the max value
        bin_end = max(y[val - number_of_gens:val]) + 1
        bins = np.arange(bin_start, bin_end, 0.5)

        ax.hist(y[val - number_of_gens:val], bins=bins)
        ax.grid(True, alpha=0.5)  # Add grid lines
        average_value = float(avg[count][-1])
        min_val = float(min(y[val - number_of_gens:val]))
        ax.axvline(average_value, color='red', linestyle='dashed',
                   linewidth=2, label='Average: ' + str(average_value))
        ax.axvline(min_val, color='green', linestyle='dashed',
                   linewidth=2, label='Min: ' + str(min_val))
        ax.legend()
        ax.set_title(c[count])

        # Append the Figure to the list
        lst.append(fig)

    return lst


def plot_line(file):
    """
    Plot line graphs for sequences from a CSV file.

    Parameters:
    - file (File): A file object containing CSV data.

    Returns:
    - lst (list): List of matplotlib figure objects.
    """
    _, y, _, amount_seq = reading_the_csv_file(file)
    lst = []

    number_of_gens = 10**5

    for count, val in enumerate(range(number_of_gens, (amount_seq + 1) * number_of_gens, number_of_gens)):
        # Create a new Figure and Axes for each plot
        fig, ax = plt.subplots()

        # Generate x values as iteration numbers
        x_values = list(range(val - number_of_gens, val))

        # Plot y values against x values
        ax.plot(x_values, y[val - number_of_gens:val])
        ax.grid(True, alpha=0.5)  # Add grid lines
        ax.set_xlabel('Iteration number')
        ax.set_ylabel('Score')
        ax.set_title(f'Sequence {count + 1}')

        # Append the Figure to the list
        lst.append(fig)

    return lst
