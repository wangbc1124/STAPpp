#!/usr/bin/env python
"""Export final-frame nodal displacements from an Abaqus ODB to CSV."""
from __future__ import print_function

import argparse
import csv
import math

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
        step_name = args.step
        if step_name is None:
            step_name = list(odb.steps.keys())[-1]
        step = odb.steps[step_name]
        frame = step.frames[args.frame]
        disp = frame.fieldOutputs["U"]

        with open(args.csv, "w") as out:
            writer = csv.writer(out)
            writer.writerow(["node", "instance", "U1", "U2", "U3", "UMag"])
            for value in disp.values:
                data = list(value.data)
                while len(data) < 3:
                    data.append(0.0)
                umag = math.sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2])
                instance = value.instance.name if value.instance is not None else ""
                writer.writerow([value.nodeLabel, instance, data[0], data[1], data[2], umag])
        print("Exported {0} nodal displacement values from step {1} frame {2}".format(
            len(disp.values), step_name, args.frame))
    finally:
        odb.close()


if __name__ == "__main__":
    main()
