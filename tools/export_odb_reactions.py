#!/usr/bin/env python
"""Export final-frame Abaqus reaction forces to CSV."""
from __future__ import print_function

import argparse
import csv
import os

from odbAccess import openOdb


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--odb", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--step", default=None)
    parser.add_argument("--frame", type=int, default=-1)
    args = parser.parse_args()

    odb = openOdb(path=args.odb, readOnly=True)
    try:
        step_name = args.step or list(odb.steps.keys())[-1]
        frame = odb.steps[step_name].frames[args.frame]
        if "RF" not in frame.fieldOutputs:
            raise RuntimeError("RF field output is not available in {}".format(args.odb))
        rf = frame.fieldOutputs["RF"]
        rows = []
        for value in rf.values:
            data = list(value.data)
            while len(data) < 3:
                data.append(0.0)
            inst = value.instance.name if value.instance is not None else ""
            node = value.instance.getNodeFromLabel(value.nodeLabel) if value.instance is not None else None
            coord = node.coordinates if node is not None else (0.0, 0.0, 0.0)
            if abs(data[0]) + abs(data[1]) + abs(data[2]) <= 0.0:
                continue
            rows.append([value.nodeLabel, inst, coord[0], coord[1], coord[2], data[0], data[1], data[2]])

        out_dir = os.path.dirname(args.csv)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.csv, "w") as out:
            writer = csv.writer(out)
            writer.writerow(["node", "instance", "x", "y", "z", "RF1", "RF2", "RF3"])
            writer.writerows(rows)
        print("Exported {0} nonzero reaction values from step {1} frame {2}".format(
            len(rows), step_name, args.frame))
    finally:
        odb.close()


if __name__ == "__main__":
    main()
