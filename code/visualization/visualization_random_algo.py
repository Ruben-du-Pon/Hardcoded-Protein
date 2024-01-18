import matplotlib.pyplot as plt
import csv
import math
import numpy as np
import os


def reading_the_csv_file(file):
    csvreader = csv.reader(file)

    rows = []
    for row in csvreader:
        rows.append(row)

    x = 1000
    c, y, avg, amount_seq = [], [], [], math.floor(len(rows)/999)

    for _ in range(amount_seq):
        for i in rows[x-1000:x]:
            y.append(int(i[-1]))
            if i[-2] not in c:
                c.append(i[-2])
        avg.append(rows[x+1])
        x += 1004

    return c, y, avg, amount_seq

def plot_random(file) -> list:
    c, y, avg, amount_seq = reading_the_csv_file(file)
    lst = []

    for count, val in enumerate(range(1000, (amount_seq+1)*1000, 1000)):
        # Create a new Figure and Axes for each plot
        fig, ax = plt.subplots()

        ax.hist(y[val-1000:val], bins=np.arange(min(y[val-1000:val]), max(y[val-1000:val]), 0.5))
        ax.grid(True, alpha=0.5)  # Add grid lines
        average_value = float(avg[count][-1])
        min_val = float(min(y[val-1000:val]))
        ax.axvline(average_value, color='red', linestyle='dashed', linewidth=2, label='Average: ' + str(average_value))
        ax.axvline(min_val, color='green', linestyle='dashed', linewidth=2, label='Min: ' + str(min_val))
        ax.legend()
        ax.set_title(c[count])

        # Append the Figure to the list
        lst.append(fig)

    return lst