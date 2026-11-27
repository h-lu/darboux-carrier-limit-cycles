"""Ordering of the specified integer bounds on the common grid domain, 5 to 350."""
import csv
from io import StringIO
from itertools import groupby


def interval_text(degrees):
    """Write maximal consecutive runs of an increasing list of integer degrees."""
    ranges = []
    for _, group in groupby(enumerate(degrees), key=lambda item: item[1] - item[0]):
        run = [degree for _, degree in group]
        ranges.append(str(run[0]) if len(run) == 1 else f"{run[0]}-{run[-1]}")
    return " ".join(ranges)


def order_table(degree_csv):
    rows = [row for row in csv.DictReader(StringIO(degree_csv))
            if 5 <= int(row["D"]) <= 350]
    if [int(row["D"]) for row in rows] != list(range(5, 351)):
        raise ValueError("comparison requires one ordered row per degree from 5 through 350")
    pairs = [("weighted", "unweighted")]
    for left in ("unweighted", "weighted", "grid"):
        for right in ("classical", "two_source", "Qin_Sun", "mechanical"):
            pairs.append((left, right))

    output = StringIO()
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["left", "right", "relation", "count", "degree_intervals"])
    for left, right in pairs:
        degrees = {key: [] for key in ("greater", "less", "equal", "unavailable")}
        for row in rows:
            if row[left] == "" or row[right] == "":
                relation = "unavailable"
            else:
                difference = int(row[left]) - int(row[right])
                relation = "equal"
                if difference > 0:
                    relation = "greater"
                elif difference < 0:
                    relation = "less"
            degrees[relation].append(int(row["D"]))
        for relation, budgets in degrees.items():
            writer.writerow([left, right, relation, len(budgets), interval_text(budgets)])
    return output.getvalue()
