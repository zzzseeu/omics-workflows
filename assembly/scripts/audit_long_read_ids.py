#!/usr/bin/env python3
"""Reject duplicate HiFi/ONT identifiers across active readsets before assembly."""
import argparse
from pathlib import Path
import sqlite3
import tempfile
import pysam
from io_utils import save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hifi',nargs='+',required=True)
    p.add_argument('--ul',nargs='+',required=True)
    p.add_argument('--out',required=True)
    args=p.parse_args()
    dest=Path(args.out);dest.parent.mkdir(parents=True,exist_ok=True)
    counts={}
    with tempfile.TemporaryDirectory(prefix='.ids-',dir=dest.parent) as temp:
        with sqlite3.connect(str(Path(temp)/'ids.sqlite')) as db:
            db.execute('PRAGMA journal_mode=OFF')
            db.execute('CREATE TABLE reads (technology TEXT, id TEXT, PRIMARY KEY(technology,id)) WITHOUT ROWID')
            for technology,paths in [('hifi',args.hifi),('ont',args.ul)]:
                count=0
                for path in paths:
                    with pysam.FastxFile(path) as reads:
                        for read in reads:
                            try:db.execute('INSERT INTO reads VALUES (?,?)',(technology,read.name))
                            except sqlite3.IntegrityError as exc:
                                raise ValueError(f'Duplicate {technology} read ID across input readsets: {read.name}') from exc
                            count+=1
                    db.commit()
                counts[technology]=count
    save(args.out,{'status':'PASS','read_counts':counts,'policy':'No silent deduplication; duplicated identifiers block assembly'})


if __name__=='__main__':main()
