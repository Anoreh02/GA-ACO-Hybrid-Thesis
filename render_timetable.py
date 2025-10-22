#!/usr/bin/env python3
"""
Render weekly timetables per room from CSV (Course,Instructor,Room,Day,Period).

This script creates one image per room (or per instructor if you choose),
showing a clean week grid (days × periods) with readable course/instructor labels.

Usage:
    python render_timetable.py timetable.csv output.png
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse

DAY_NAMES = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]

def parse_args():
    p = argparse.ArgumentParser(description="Render per-room weekly timetable(s).")
    p.add_argument("csv", help="Input timetable.csv file")
    p.add_argument("out", help="Output file prefix, e.g., week.png → week_room-0.png etc.")
    p.add_argument("--days", type=int, default=5, help="Number of days (default=5)")
    p.add_argument("--periods", type=int, default=6, help="Number of periods per day (default=6)")
    p.add_argument("--group-by", choices=["room", "instructor"], default="room",
                   help="Generate timetables grouped by room (default) or instructor.")
    p.add_argument("--dpi", type=int, default=150, help="Image DPI")
    return p.parse_args()

def build_grid(df, D, P):
    """Build a P×D grid for one room/instructor."""
    grid = [["" for _ in range(D)] for __ in range(P)]
    for _, row in df.iterrows():
        d = int(row["Day"]) - 1
        p = int(row["Period"]) - 1
        course = str(row["Course"])
        instr = str(row["Instructor"])
        room = str(row["Room"])
        label = f"Course {course}\nInst {instr}"
        if grid[p][d]:
            grid[p][d] += f"\n—\n{label}"
        else:
            grid[p][d] = label
    return grid

def draw_table(grid, title, D, P, filename, dpi):
    """Draw and save a weekly table image."""
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axis("off")

    day_labels = DAY_NAMES[:D]
    period_labels = [f"P{p+1}" for p in range(P)]

    # Build text table (rows × cols)
    nrows, ncols = P + 1, D + 1
    table_text = [["" for _ in range(ncols)] for __ in range(nrows)]
    for r in range(nrows):
        for c in range(ncols):
            if r == 0 and c == 0:
                table_text[r][c] = ""
            elif r == 0:
                table_text[r][c] = day_labels[c - 1]
            elif c == 0:
                table_text[r][c] = period_labels[r - 1]
            else:
                table_text[r][c] = grid[r - 1][c - 1]

    table = ax.table(cellText=table_text, loc="center", cellLoc="center", colLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.4)

    # Styling
    for (r, c), cell in table.get_celld().items():
        if r == 0 or c == 0:
            cell.set_facecolor("#e0e0e0")
            cell.set_fontsize(9)
            cell.set_text_props(weight="bold")
        else:
            cell.set_facecolor("#f9f9f9")

    ax.set_title(title, fontsize=13, pad=10)
    fig.tight_layout()
    fig.savefig(filename, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved {filename}")

def main():
    args = parse_args()
    df = pd.read_csv(args.csv)
    required = ["Course", "Instructor", "Room", "Day", "Period"]
    if not all(col in df.columns for col in required):
        raise ValueError(f"CSV must have columns: {required}")

    D, P = args.days, args.periods

    if args.group_by == "room":
        groups = sorted(df["Room"].unique(), key=lambda x: str(x))
        keycol, label = "Room", "Room"
        suffix = "room"
    else:
        groups = sorted(df["Instructor"].unique(), key=lambda x: str(x))
        keycol, label = "Instructor", "Instructor"
        suffix = "inst"

    base = args.out.rsplit(".", 1)[0]
    for g in groups:
        sub = df[df[keycol] == g]
        grid = build_grid(sub, D, P)
        filename = f"{base}_{suffix}-{g}.png"
        title = f"{label} {g}"
        draw_table(grid, title, D, P, filename, args.dpi)

if __name__ == "__main__":
    main()
