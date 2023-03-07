#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import List, Dict, Tuple
import argparse
import logging


def arguments() -> Tuple[str, str, str]:
    # initialize parser
    msg = "Adding options"
    parser = argparse.ArgumentParser(description=msg)
    # add some options
    parser.add_argument(
        "-sdf",
        "--sdf",
        type=str,
        default="",
        help="input mol.sdf",
    )
    parser.add_argument(
        "-sdf_folder",
        "--sdf_folder",
        type=str,
        default=".",
        help="the folder to store multiple sdf files",
    )
    parser.add_argument(
        "-build_fingerprint_db",
        "--build_fingerprint_db",
        action="store_true",
        help="build fingerprint catalog or not",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        default="results",
        help="output folder to store fp_catalogue",
    )
    parser.add_argument(
        "-load_fp_catalogue",
        "--load_fp_catalogue",
        type=str,
        default=".",
        help="output folder to check if filter workflow already done or not e.g. /fsx/bwamem/SAMPLE_filtered.bam",
    )
    parser.add_argument(
        "-threads",
        "--threads",
        type=int,
        default=1,
        help="threads to generate fragments catalog",
    )
    # read in arguments and print
    args = parser.parse_args()

    sdf = args.sdf
    sdf_folder = args.sdf_folder
    build_fingerprint_db = args.build_fingerprint_db
    outdir = args.outdir
    threads = args.threads

    return sdf, sdf_folder, build_fingerprint_db, outdir, threads


def setup_logger(name: str = "Default", logfile: str = None) -> logging.Logger:
    logging.getLogger().setLevel(1)
    logger_string_format = "[%(name)s - %(asctime)s] %(message)s"
    logging.basicConfig(format=logger_string_format)
    logger = logging.getLogger(name)

    if logfile is not None:
        hdlr = logging.FileHandler(logfile)
        formatter = logging.Formatter(logger_string_format)
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr)
    return logger
