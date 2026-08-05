#!/usr/bin/env python3
"""Prepare compact, deterministic inputs for dynamic Swiggum cluster snapping."""

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path


STAR_OUTPUT = "sn_history_july9_imf_sample0_snap_ready_minus100_to_plus20myr.csv.gz"
CLUSTER_OUTPUT = "cluster_births_july9_imf_sample0.csv.gz"
AUDIT_OUTPUT = "cluster_births_july9_imf_sample0.audit.json"


def _open_text(path, mode):
    if str(path).endswith(".gz"):
        if "w" in mode:
            raw = open(path, "wb")
            zipped = gzip.GzipFile(fileobj=raw, mode="wb", mtime=0)
            return io.TextIOWrapper(zipped, encoding="utf-8", newline="")
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def _read_csv(path):
    with _open_text(path, "rt") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV file has no header: {path}")
        return list(reader.fieldnames), list(reader)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _one_float(rows, column, cluster_name):
    values = {float(row[column]) for row in rows}
    if len(values) != 1:
        raise ValueError(
            f"Ambiguous {column} values for cluster {cluster_name!r}: {sorted(values)}"
        )
    value = values.pop()
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {column} for cluster {cluster_name!r}")
    return value


def _interpolate_birth(rows, birth_time, cluster_name):
    by_time = {}
    for row in rows:
        time = float(row["time_myr"])
        position = tuple(
            float(row[key])
            for key in ("X_cluster_pc", "Y_cluster_pc", "Z_cluster_pc")
        )
        if time in by_time and by_time[time] != position:
            raise ValueError(
                f"Ambiguous trajectory positions for {cluster_name!r} at {time} Myr"
            )
        by_time[time] = position
    times = sorted(by_time)
    if not times or birth_time < times[0] or birth_time > times[-1]:
        raise ValueError(
            f"Birth time {birth_time} Myr for {cluster_name!r} is outside "
            f"trajectory range [{times[0]}, {times[-1]}]"
        )
    if birth_time in by_time:
        return by_time[birth_time], {
            "lower_time_myr": birth_time,
            "upper_time_myr": birth_time,
            "fraction": 0.0,
            "exact_sample": True,
        }
    upper_index = next(i for i, time in enumerate(times) if time > birth_time)
    lower_time, upper_time = times[upper_index - 1], times[upper_index]
    fraction = (birth_time - lower_time) / (upper_time - lower_time)
    lower, upper = by_time[lower_time], by_time[upper_time]
    position = tuple(
        lower[i] + fraction * (upper[i] - lower[i]) for i in range(3)
    )
    return position, {
        "lower_time_myr": lower_time,
        "upper_time_myr": upper_time,
        "fraction": fraction,
        "exact_sample": False,
    }


def prepare(stellar_path, cluster_path, output_dir, history_start_time_myr,
            stellar_output=STAR_OUTPUT, cluster_output=CLUSTER_OUTPUT,
            audit_output=AUDIT_OUTPUT):
    stellar_path, cluster_path = Path(stellar_path), Path(cluster_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    star_header, star_rows = _read_csv(stellar_path)
    cluster_header, cluster_rows = _read_csv(cluster_path)

    required_star = {"star_id", "star_cluster_name"}
    required_cluster = {
        "cluster_name", "time_myr", "X_cluster_pc", "Y_cluster_pc",
        "Z_cluster_pc", "cluster_age_myr", "cluster_radius_pc",
    }
    missing_star = sorted(required_star - set(star_header))
    missing_cluster = sorted(required_cluster - set(cluster_header))
    if missing_star or missing_cluster:
        raise ValueError(
            f"Missing columns: stellar={missing_star}, cluster={missing_cluster}"
        )
    if "snap_cluster_id" in star_header:
        raise ValueError("Input stellar history already has snap_cluster_id")

    cluster_rows_by_name = defaultdict(list)
    for row in cluster_rows:
        cluster_rows_by_name[row["cluster_name"]].append(row)
    star_names_by_id = defaultdict(set)
    members_by_name = defaultdict(set)
    for row in star_rows:
        name, star_id = row["star_cluster_name"], row["star_id"]
        star_names_by_id[star_id].add(name)
        members_by_name[name].add(star_id)
    ambiguous_stars = {
        star_id: sorted(names) for star_id, names in star_names_by_id.items()
        if len(names) != 1
    }
    referenced_names = sorted(members_by_name)
    cluster_names = set(cluster_rows_by_name)
    unmatched_stars = sorted(set(referenced_names) - cluster_names)
    unmatched_clusters = sorted(cluster_names - set(referenced_names))
    unmatched_star_ids = sorted(
        {
            (star_id, next(iter(names)))
            for star_id, names in star_names_by_id.items()
            if len(names) == 1 and next(iter(names)) in unmatched_stars
        }
    )

    audit = {
        "inputs": {
            "stellar_history": str(stellar_path.resolve()),
            "cluster_history": str(cluster_path.resolve()),
        },
        "history_start_time_myr": history_start_time_myr,
        "counts": {
            "stellar_rows": len(star_rows),
            "unique_stars": len(star_names_by_id),
            "referenced_clusters": len(referenced_names),
            "cluster_history_clusters": len(cluster_names),
            "matched_clusters": len(set(referenced_names) & cluster_names),
            "matched_stars": sum(
                len(members_by_name[name])
                for name in set(referenced_names) & cluster_names
            ),
        },
        "matched_cluster_names": sorted(set(referenced_names) & cluster_names),
        "unmatched_stars": [
            {"star_id": star_id, "star_cluster_name": name}
            for star_id, name in unmatched_star_ids
        ],
        "unmatched_star_clusters": unmatched_stars,
        "unmatched_cluster_history_clusters": unmatched_clusters,
        "ambiguous_star_cluster_assignments": ambiguous_stars,
        "interpolation_checks": [],
        "sha256": {
            "stellar_input": _sha256(stellar_path),
            "cluster_input": _sha256(cluster_path),
        },
    }
    audit_path = output_dir / audit_output
    if ambiguous_stars or unmatched_stars:
        audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
        raise ValueError(
            "Exact star_cluster_name == cluster_name join failed; see "
            f"{audit_path}"
        )

    ids = {name: index for index, name in enumerate(referenced_names)}
    birth_rows = []
    for name in referenced_names:
        rows = cluster_rows_by_name[name]
        age = _one_float(rows, "cluster_age_myr", name)
        radius = _one_float(rows, "cluster_radius_pc", name)
        birth_time = -age
        if not birth_time > history_start_time_myr:
            raise ValueError(
                f"Cluster {name!r} birth {birth_time} Myr is not later than "
                f"history start {history_start_time_myr} Myr"
            )
        position, check = _interpolate_birth(rows, birth_time, name)
        audit["interpolation_checks"].append({
            "snap_cluster_id": ids[name], "cluster_name": name,
            "birth_time_myr": birth_time, **check,
        })
        birth_rows.append({
            "snap_cluster_id": ids[name], "cluster_name": name,
            "birth_time_myr": format(birth_time, ".17g"),
            "birth_x_pc": format(position[0], ".17g"),
            "birth_y_pc": format(position[1], ".17g"),
            "birth_z_pc": format(position[2], ".17g"),
            "cluster_radius_pc": format(radius, ".17g"),
            "number_of_sampled_stars": len(members_by_name[name]),
        })

    star_output_path = output_dir / stellar_output
    with _open_text(star_output_path, "wt") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=star_header + ["snap_cluster_id"],
            lineterminator="\n",
        )
        writer.writeheader()
        for row in star_rows:
            writer.writerow(
                {**row, "snap_cluster_id": ids[row["star_cluster_name"]]}
            )
    cluster_output_path = output_dir / cluster_output
    birth_header = [
        "snap_cluster_id", "cluster_name", "birth_time_myr", "birth_x_pc",
        "birth_y_pc", "birth_z_pc", "cluster_radius_pc",
        "number_of_sampled_stars",
    ]
    with _open_text(cluster_output_path, "wt") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=birth_header, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(birth_rows)

    audit["counts"]["output_clusters"] = len(birth_rows)
    audit["counts"]["output_stellar_rows"] = len(star_rows)
    audit["sha256"]["stellar_output"] = _sha256(star_output_path)
    audit["sha256"]["cluster_birth_output"] = _sha256(cluster_output_path)
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    return star_output_path, cluster_output_path, audit_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stellar_history",
        nargs="?",
        default="sn_history_july9_imf_sample0_minus100_to_plus20myr.csv.gz",
    )
    parser.add_argument(
        "cluster_history",
        nargs="?",
        default="cluster_history_july9_imf_sample0_minus100_to_plus20myr.csv.gz",
    )
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--stellar-output", default=STAR_OUTPUT)
    parser.add_argument("--cluster-output", default=CLUSTER_OUTPUT)
    parser.add_argument("--audit-output", default=AUDIT_OUTPUT)
    parser.add_argument("--history-start-time-myr", type=float, default=-100.0)
    args = parser.parse_args()
    outputs = prepare(args.stellar_history, args.cluster_history,
                      args.output_dir, args.history_start_time_myr,
                      args.stellar_output, args.cluster_output,
                      args.audit_output)
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
