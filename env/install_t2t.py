#!/usr/bin/env python3
"""Install pinned post-assembly tools additively inside the unified prefix.

Run on a compute node with that prefix's Python. No conda solve or base-package
upgrade is performed. Downloads and source trees live outside the workflow repo.
"""
import argparse
import hashlib
import json
import subprocess
import sys
import tarfile
import urllib.request
import zipfile
from pathlib import Path

SOURCES = {
    "quartet": ("aaranyue/quarTeT", "a44d8c87e3f4bcefaf43dc4134e88b6a21890265", "1.3.1"),
    "inspector": ("Maggi-Chen/Inspector", "2d813d96eb38adec9b8fb25388a3f4f1a2c63573", "1.3"),
}
NP2_URL = "https://conda.anaconda.org/bioconda/linux-64/nextpolish2-0.2.2-h74ec884_0.conda"
NP2_SHA = "87a1ff84874c3f036f2440d919d3f75ce0169c9a49c3820751a21ef1bceec98e"
TIDK_URL = "https://conda.anaconda.org/bioconda/linux-64/tidk-0.3.1-hef77e2b_0.conda"
TIDK_SHA = "6f4ef177e317193c3ad1ae132534b9b8947440a849c275f356e0b6fe4abf6570"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument(
        "--download-only",
        action="store_true",
        help="Prefetch on a network-enabled host, without changing the prefix",
    )
    args = p.parse_args()
    if args.download_only:
        args.cache.mkdir(parents=True, exist_ok=True)
        items = [
            (f"https://codeload.github.com/{repo}/tar.gz/{commit}", name + ".tar.gz")
            for name, (repo, commit, version) in SOURCES.items()
        ]
        items.extend([(NP2_URL, "nextPolish2.conda"), (TIDK_URL, "tidk.conda")])
        for url, name in items:
            if not (args.cache / name).exists():
                urllib.request.urlretrieve(url, args.cache / name)
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "download",
                "--no-deps",
                "--dest",
                str(args.cache),
                "statsmodels==0.15.0",
                "patsy==1.0.2",
            ],
            check=True,
        )
        print("T2T_DOWNLOADS=READY")
        return
    prefix = Path(sys.prefix).resolve()
    root = prefix / "opt/T2T"
    if (root / "SOURCE_LOCK.json").exists():
        raise SystemExit("T2T installation already locked; verify it instead of overwriting")
    before = sha(prefix / "conda-meta/history")
    root.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)
    (root / "bin").mkdir(exist_ok=True)
    archives = {}

    def download(url, name):
        path = args.cache / name
        if not path.exists():
            urllib.request.urlretrieve(url, path)
        archives[name] = {"url": url, "sha256": sha(path)}
        return path

    for name, (repo, commit, version) in SOURCES.items():
        archive = download(f"https://codeload.github.com/{repo}/tar.gz/{commit}", name + ".tar.gz")
        destination = root / name
        if not destination.exists():
            destination.mkdir()
            with tarfile.open(archive) as tar:
                for member in tar.getmembers():
                    member.name = "/".join(member.name.split("/")[1:])
                    if member.name:
                        tar.extract(member, destination, filter="data")

    def conda_binary(url, name, expected):
        archive = download(url, name + ".conda")
        if sha(archive) != expected:
            raise ValueError(name + " package checksum mismatch")
        with zipfile.ZipFile(archive) as z:
            member = next(n for n in z.namelist() if n.startswith("pkg-"))
            compressed = args.cache / (name + ".pkg.tar.zst")
            compressed.write_bytes(z.read(member))
        # Explicit cleanup below drains the streaming extractor.
        # pylint: disable-next=consider-using-with
        process = subprocess.Popen(
            [str(prefix / "bin/zstd"), "-d", "-c", str(compressed)], stdout=subprocess.PIPE
        )
        try:
            with tarfile.open(fileobj=process.stdout, mode="r|") as tar:
                for member in tar:
                    if member.name == "bin/" + name:
                        tar.extract(member, root, filter="data")
            if process.wait():
                raise RuntimeError("zstd extraction failed")
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
                process.wait()

    conda_binary(NP2_URL, "nextPolish2", NP2_SHA)
    conda_binary(TIDK_URL, "tidk", TIDK_SHA)
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-index",
            "--no-deps",
            "--find-links",
            str(args.cache),
            "--target",
            str(root / "python"),
            "statsmodels==0.15.0",
            "patsy==1.0.2",
        ],
        check=True,
    )
    wrappers = {
        "quartet": (
            f'import sys\nsys.path.insert(0, {str(root/"quartet")!r})\n'
            "from quartet.main import main\nmain()\n"
        ),
        "inspector.py": (
            f'import runpy, sys\nsys.path.insert(0, {str(root/"python")!r})\n'
            f'sys.path.insert(0, {str(root/"inspector")!r})\n'
            f'runpy.run_path({str(root/"inspector/inspector.py")!r}, run_name="__main__")\n'
        ),
    }
    (root / "libexec").mkdir(exist_ok=True)
    for name in ("tidk", "nextPolish2"):
        native = root / "libexec" / name
        (root / "bin" / name).replace(native)
        native.chmod(0o755)
        library = str(prefix / "lib")
        wrappers[name] = (
            f"import os, sys\nenv=dict(os.environ)\n"
            f'env["LD_LIBRARY_PATH"]={library!r}+os.pathsep+env.get("LD_LIBRARY_PATH", "")\n'
            f"os.execve({str(native)!r}, [{str(native)!r}, *sys.argv[1:]], env)\n"
        )
    for name, body in wrappers.items():
        (root / "bin" / name).write_text(f"#!{prefix}/bin/python\n" + body)
    for name in ["quartet", "inspector.py", "tidk", "nextPolish2"]:
        executable = root / "bin" / name
        executable.chmod(0o755)
        link = prefix / "bin" / name
        if link.is_symlink() and link.resolve() == executable:
            continue
        if link.exists() or link.is_symlink():
            if name == "tidk" and link.resolve().is_relative_to(prefix):
                version = subprocess.run([str(link), "--version"], capture_output=True, check=False)
                if version.returncode == 0 and version.stdout.decode().strip() == "tidk 0.3.1":
                    continue
            raise ValueError(f"Refusing to overwrite {link}")
        link.symlink_to(executable)
    for name in ["quartet", "inspector.py", "tidk", "nextPolish2"]:
        result = subprocess.run([str(root / "bin" / name), "--help"], capture_output=True, check=False)
        if result.returncode:
            raise RuntimeError(f'{name} help failed: {result.stderr.decode(errors="replace")}')
    if sha(prefix / "conda-meta/history") != before:
        raise ValueError("Base environment history changed")
    lock = {
        "tools": {
            **{k: {"repository": v[0], "commit": v[1], "version": v[2]} for k, v in SOURCES.items()},
            "tidk": {"version": "0.3.1"},
            "nextPolish2": {"version": "0.2.2"},
        },
        "archives": archives,
        "base_history_sha256": before,
        "files": {
            str(f.relative_to(root)): sha(f)
            for f in sorted(root.rglob("*"))
            if f.is_file() and "__pycache__" not in f.parts
        },
    }
    (root / "SOURCE_LOCK.json").write_text(json.dumps(lock, sort_keys=True, indent=2) + "\n")
    print("T2T_ADDITIVE_INSTALL=PASS", flush=True)


if __name__ == "__main__":
    main()
