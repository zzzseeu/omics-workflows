"""Lossless read preparation and verified reuse of historical prepared inputs."""

import csv
import fcntl
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from assembly.lib.io import checksum, execute, load, save, signature, source_check

POLICY = "lossless-prepared-fastq-reuse-v1"


def metadata(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size > 1024 * 1024:
        raise ValueError(f"Missing or oversized prepared-input metadata: {path}")
    content = path.read_bytes()
    return json.loads(content), hashlib.sha256(content).hexdigest()


def describe_prepared_input(spec, manifest, cache_dir, run_dir, repository):
    """No reads, decompression or large-file hashing during preflight."""
    sid, rid = spec["sample_id"], spec["readset_id"]
    if spec["technology"] not in ("hifi", "ont"):
        raise ValueError("Prepared FASTQ reuse is supported only for HiFi/ONT")
    manifest, cache_dir, run_dir, repository = map(
        lambda p: Path(p).resolve(), (manifest, cache_dir, run_dir, repository)
    )
    for path in (manifest, cache_dir):
        if (
            path == run_dir
            or path.is_relative_to(run_dir)
            or path == repository
            or path.is_relative_to(repository)
        ):
            raise ValueError(
                "Prepared-input metadata and cache must be external to the new run and workflow repository"
            )
    if cache_dir.exists() and not cache_dir.is_dir():
        raise ValueError("Prepared-input cache is not a directory")
    old, manifest_sha = metadata(manifest)
    for key in ("sample_id", "readset_id", "library_id"):
        if old.get(key) != spec[key]:
            raise ValueError("Prepared-input identity mismatch: " + key)
    if old.get("technology", spec["technology"]) != spec["technology"]:
        raise ValueError("Prepared-input technology mismatch")
    if Path(old.get("source", "")).resolve() != Path(spec["file1"]).resolve():
        raise ValueError("Prepared-input raw source path mismatch")
    if old.get("source_signature") != spec["source_signature"]:
        raise ValueError("Prepared-input raw source changed after preparation")
    if old.get("filtering") != "No quality trimming or deduplication performed":
        raise ValueError("Prepared-input filtering policy is incompatible")
    for key, size in [("source_md5", 32), ("source_sha256", 64)]:
        if not re.fullmatch("[0-9a-f]{" + str(size) + "}", old.get(key, "")):
            raise ValueError("Prepared-input source checksum missing: " + key)
    if spec.get("expected_md5") and old["source_md5"] != spec["expected_md5"]:
        raise ValueError("Prepared-input raw source MD5 mismatch")
    if old.get("source_format", spec["format"]) != spec["format"]:
        raise ValueError("Prepared-input source format mismatch")
    conversion = old.get("conversion", "")
    if not old.get("reuse") and not conversion.startswith(
        "samtools fastq" if spec["format"] == "bam" else "FASTQ byte copy"
    ):
        raise ValueError("Prepared-input conversion policy is incompatible")
    parent = manifest.parent
    if parent.name != "01_inputs" or parent.parent.name != sid or parent.parent.parent.name != "samples":
        raise ValueError(
            "Prepared-input manifest must belong to a numbered samples/<sample>/01_inputs directory"
        )
    expected = parent / (rid + ".fastq.gz")
    declared = Path(old.get("prepared_fastq", ""))
    if declared.is_absolute():
        if declared.resolve() != expected.resolve():
            raise ValueError("Prepared FASTQ path disagrees with its manifest")
    elif declared.as_posix() != f"samples/{sid}/01_inputs/{rid}.fastq.gz":
        raise ValueError("Unsupported relative prepared FASTQ path")
    fastq = expected.resolve()
    if (
        not fastq.is_file()
        or not fastq.stat().st_size
        or fastq.is_relative_to(run_dir)
        or fastq.is_relative_to(repository)
    ):
        raise ValueError("Prepared FASTQ is missing, empty or not external")
    qc_path = (
        Path(old["reuse"]["qc_manifest"])
        if old.get("reuse")
        else parent.parent / "02_read_qc" / (rid + ".json")
    )
    qc, qc_sha = metadata(qc_path)
    if old.get("reuse") and old["reuse"].get("qc_manifest_sha256") != qc_sha:
        raise ValueError("Inherited prepared-input QC checksum mismatch")
    if qc.get("status") != "PASS" or qc.get("num_reads", 0) < 1 or qc.get("num_bases", 0) < 1:
        raise ValueError("Prepared-input QC is not complete/PASS")
    recorded = qc.get("statistics", {}).get("file", "")
    expected_relative = f"samples/{sid}/01_inputs/{rid}.fastq.gz"
    if (
        not recorded
        or (Path(recorded).is_absolute() and Path(recorded).resolve() != fastq)
        or (not Path(recorded).is_absolute() and recorded != expected_relative)
    ):
        raise ValueError("Prepared-input QC refers to a different FASTQ")
    observed = signature(fastq)
    producer_sha = old.get("prepared_sha256")
    if producer_sha is not None:
        if not re.fullmatch("[0-9a-f]{64}", producer_sha) or not isinstance(
            old.get("prepared_signature"), dict
        ):
            raise ValueError("Prepared FASTQ producer checksum/signature record is invalid")
        if old["prepared_signature"].get("size") != observed["size"]:
            raise ValueError("Prepared FASTQ size differs from its producer record")
    identity = {
        "policy": POLICY,
        "sample_id": sid,
        "readset_id": rid,
        "library_id": spec["library_id"],
        "technology": spec["technology"],
        "prepared_fastq": str(fastq),
        "manifest": str(manifest),
        "manifest_sha256": manifest_sha,
        "qc": str(qc_path),
        "qc_sha256": qc_sha,
        "source_sha256": old["source_sha256"],
    }
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    return {
        "identity": identity,
        "prepared_signature": observed,
        "manifest_signature": signature(manifest),
        "qc_signature": signature(qc_path),
        "producer_sha256": producer_sha,
        "expected_num_reads": qc["num_reads"],
        "expected_num_bases": qc["num_bases"],
        "producer_signature_matches": old.get("prepared_signature") == observed,
        "cache_file": str(cache_dir / (key + ".json")),
    }


def prepare_long_reads(args):
    spec = load(args.spec)
    source = Path(spec["file1"])
    before = signature(source)
    if before != spec["source_signature"]:
        raise ValueError("Source changed after preflight; restart the launcher to create a new plan")
    md5, sha = hashlib.md5(), hashlib.sha256()
    with source.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5.update(block)
            sha.update(block)
    if spec["expected_md5"] and md5.hexdigest() != spec["expected_md5"]:
        raise ValueError(f"MD5 mismatch: {source}; observed {md5.hexdigest()}")
    output = Path(args.fastq)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".partial.fastq.gz")
    try:
        if spec["format"] == "bam":
            execute(["samtools", "quickcheck", "-u", "-v", source])
            paired = (
                execute(["samtools", "view", "-c", "-f", "1", source], stdout=subprocess.PIPE)
                .stdout.decode()
                .strip()
            )
            if int(paired):
                raise ValueError("HiFi BAM contains paired reads; only unpaired HiFi BAM is supported")
            # Unpaired primary reads go to category 0. Native compression also works
            # with a one-core budget without a second Python compressor process.
            execute(
                [
                    "samtools",
                    "fastq",
                    "-n",
                    "-@",
                    str(max(0, args.threads - 1)),
                    "-0",
                    temporary,
                    "-s",
                    os.devnull,
                    source,
                ],
                stdout=subprocess.DEVNULL,
            )
            conversion = "samtools fastq; unpaired primary reads; original IDs retained"
        else:
            with source.open("rb") as probe:
                compressed = probe.read(2) == b"\x1f\x8b"
            if compressed:
                temporary.symlink_to(source.resolve())
            else:
                with source.open("rb") as src, gzip.open(temporary, "wb", compresslevel=3) as target:
                    shutil.copyfileobj(src, target, 8 * 1024 * 1024)
            conversion = "Immutable FASTQ link (gzip) or lossless gzip compression (plain text)"
        if signature(source) != before:
            raise ValueError("Source changed while preparing reads")
        # Compressed FASTQ byte copies already have exactly the source hash.
        # Converted/compressed outputs get a producer hash for future reuse.
        if spec["format"] == "fastq" and compressed:
            prepared_sha = sha.hexdigest()
        else:
            prepared = hashlib.sha256()
            with temporary.open("rb") as handle:
                for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                    prepared.update(block)
            prepared_sha = prepared.hexdigest()
        temporary.replace(output)
        save(
            args.manifest,
            {
                "sample_id": spec["sample_id"],
                "readset_id": spec["readset_id"],
                "library_id": spec["library_id"],
                "source": str(source),
                "source_signature": before,
                "source_md5": md5.hexdigest(),
                "source_sha256": sha.hexdigest(),
                "provided_md5_status": "PASS" if spec["expected_md5"] else "NOT_PROVIDED",
                "preparation_format_version": 2,
                "technology": spec["technology"],
                "source_format": spec["format"],
                "prepared_fastq": str(output.absolute()),
                "conversion": conversion,
                "prepared_sha256": prepared_sha,
                "prepared_signature": signature(output),
                "filtering": "No quality trimming or deduplication performed",
            },
        )
    finally:
        temporary.unlink(missing_ok=True)


def unchanged(spec):
    source_check(spec)
    r = spec["prepared_reuse"]
    identity = r["identity"]
    for key, sigkey, hashkey in [
        ("manifest", "manifest_signature", "manifest_sha256"),
        ("qc", "qc_signature", "qc_sha256"),
    ]:
        path = identity[key]
        if signature(path) != r[sigkey] or checksum(path) != identity[hashkey]:
            raise ValueError("Prepared-input metadata changed after preflight: " + key)
    if signature(identity["prepared_fastq"]) != r["prepared_signature"]:
        raise ValueError("Prepared FASTQ changed after preflight")


def scan_legacy(fastq, seqkit, expected_reads, expected_bases, folder):
    """Hash compressed bytes and stream once through native seqkit FASTQ parsing."""
    sha = hashlib.sha256()
    with tempfile.TemporaryFile(mode="w+b", dir=folder) as stats:
        command = [seqkit, "stats", "-T", "-j", "1", "-"]
        print("VERIFY_LEGACY_PREPARED_FASTQ", fastq, flush=True)
        # Explicit process cleanup below owns the pipe/signal lifetime.
        # pylint: disable-next=consider-using-with
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=stats)
        try:
            with Path(fastq).open("rb") as src:
                for block in iter(lambda: src.read(8 * 1024 * 1024), b""):
                    sha.update(block)
                    process.stdin.write(block)
            process.stdin.close()
            if process.wait():
                raise ValueError("Prepared FASTQ parsing failed")
        finally:
            if process.stdin and not process.stdin.closed:
                process.stdin.close()
            if process.poll() is None:
                process.terminate()
                process.wait()
        stats.seek(0)
        rows = list(csv.DictReader(stats.read().decode().splitlines(), delimiter="\t"))
        if len(rows) != 1 or rows[0]["format"].upper() != "FASTQ":
            raise ValueError("Prepared input is not FASTQ")
        reads = int(rows[0]["num_seqs"].replace(",", ""))
        bases = int(rows[0]["sum_len"].replace(",", ""))
        if (reads, bases) != (expected_reads, expected_bases):
            raise ValueError("Prepared FASTQ does not match prior QC counts")
    return sha.hexdigest()


def attest(spec):
    unchanged(spec)
    r = spec["prepared_reuse"]
    identity = r["identity"]
    cache = Path(r["cache_file"])
    cache.parent.mkdir(parents=True, exist_ok=True)
    with cache.with_suffix(".lock").open("a", encoding="utf-8") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        unchanged(spec)
        hit = cache.exists()
        if hit:
            if cache.is_symlink():
                raise ValueError("Prepared-input cache record must not be a symlink")
            record = load(cache)
            if (
                record.get("identity") != identity
                or record.get("status") != "PASS"
                or record.get("prepared_signature") != r["prepared_signature"]
            ):
                raise ValueError("Cached prepared FASTQ changed or cache provenance is invalid")
            if not re.fullmatch("[0-9a-f]{64}", record.get("prepared_sha256", "")):
                raise ValueError("Prepared-input cache checksum is invalid")
            if r["producer_sha256"] and record["prepared_sha256"] != r["producer_sha256"]:
                raise ValueError("Prepared-input cache checksum disagrees with producer")
        else:
            if r["producer_sha256"] and r["producer_signature_matches"]:
                sha = r["producer_sha256"]
                method = "PRODUCER_SHA256_SIGNATURE_AND_PASS_QC"
            else:
                sha = scan_legacy(
                    identity["prepared_fastq"],
                    spec["seqkit"]["path"],
                    r["expected_num_reads"],
                    r["expected_num_bases"],
                    cache.parent,
                )
                if r["producer_sha256"] and sha != r["producer_sha256"]:
                    raise ValueError("Prepared FASTQ producer SHA256 mismatch")
                method = (
                    "PRODUCER_SHA256_REVERIFIED_AFTER_TIMESTAMP_CHANGE"
                    if r["producer_sha256"]
                    else "LEGACY_ONE_PASS_COMPRESSED_SHA256_AND_FASTQ_QC_COUNTS"
                )
            unchanged(spec)
            record = {
                "status": "PASS",
                "identity": identity,
                "prepared_signature": r["prepared_signature"],
                "prepared_sha256": sha,
                "validation": method,
                "num_reads": r["expected_num_reads"],
                "num_bases": r["expected_num_bases"],
                "created_utc": datetime.now(timezone.utc).isoformat(),
            }
            save(cache, record)
    unchanged(spec)
    return record, hit, checksum(cache)


def reuse(args):
    spec = load(args.spec)
    record, hit, cache_sha = attest(spec)
    r = spec["prepared_reuse"]
    source = Path(r["identity"]["prepared_fastq"])
    output = Path(args.fastq)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.resolve() == source:
        raise ValueError("Reuse destination aliases the original FASTQ")
    old = load(r["identity"]["manifest"])
    temporary = output.with_name(output.name + ".reuse.partial")
    try:
        if temporary.exists() or temporary.is_symlink():
            temporary.unlink()
        temporary.symlink_to(source)
        unchanged(spec)
        temporary.replace(output)
        # Preserve the original raw-source provenance; never replace the raw
        # source hash with the compressed FASTQ hash.
        save(
            args.manifest,
            {
                **{
                    k: old[k]
                    for k in (
                        "sample_id",
                        "readset_id",
                        "library_id",
                        "source",
                        "source_signature",
                        "source_md5",
                        "source_sha256",
                        "provided_md5_status",
                        "conversion",
                        "filtering",
                    )
                },
                "preparation_format_version": 2,
                "technology": spec["technology"],
                "source_format": spec["format"],
                "prepared_fastq": str(output.absolute()),
                "prepared_sha256": record["prepared_sha256"],
                "prepared_signature": record["prepared_signature"],
                "reuse": {
                    "policy": POLICY,
                    "original_manifest": r["identity"]["manifest"],
                    "original_manifest_sha256": r["identity"]["manifest_sha256"],
                    "qc_manifest": r["identity"]["qc"],
                    "qc_manifest_sha256": r["identity"]["qc_sha256"],
                    "prepared_target": str(source),
                    "cache_file": r["cache_file"],
                    "cache_sha256": cache_sha,
                    "validation": record["validation"],
                    "cache_hit": hit,
                    "storage": "SYMLINK_NO_COPY",
                },
            },
        )
    finally:
        if temporary.exists() or temporary.is_symlink():
            temporary.unlink()
    print(
        "PREPARED_INPUT_REUSED",
        spec["sample_id"],
        spec["readset_id"],
        "CACHE_HIT" if hit else "ATTESTED",
        flush=True,
    )
