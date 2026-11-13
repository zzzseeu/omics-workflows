"""Atomic publication, integrity checks and external process execution."""

import hashlib
import json
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path


def load(path):
    """从UTF-8 JSON文件读取并返回Python对象；格式或文件错误直接抛出。"""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def save(path, value):
    """将对象按固定格式写为JSON；自动创建父目录，先写临时文件再替换目标。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)


def signature(path):
    """返回文件大小、时间和inode组成的轻量签名，用于判断是否需要重新计算完整校验值。"""
    stat = Path(path).stat()
    return {
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
        "inode": stat.st_ino,
    }


def checksum(path, algorithm="sha256"):
    """按指定哈希算法分块读取文件，返回十六进制校验值；默认SHA-256。"""
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def execute(command, **kwargs):
    """以参数列表执行外部工具并记录命令，返回CompletedProcess；非零退出码抛出异常。"""
    print("COMMAND", json.dumps([str(x) for x in command]), flush=True)
    return subprocess.run([str(x) for x in command], check=True, **kwargs)


def source_check(spec):
    """对照预检记录核验单端或双端输入的轻量签名，防止运行期间替换原始数据。"""
    for i in (1, 2):
        key = "source_signature" if i == 1 else "source_signature_2"
        if spec.get(f"file{i}") and signature(spec[f"file{i}"]) != spec[key]:
            raise ValueError(f'Source changed after preflight: {spec[f"file{i}"]}')


def publish(stage, destination, scratch):
    """把暂存结果发布到目标目录：先备份受管旧目录，再移动新目录，拒绝覆盖非受管文件。"""
    destination = Path(destination)
    if destination.exists():
        # Snakemake pre-creates every output parent, including nested empty
        # directories. Only actual files/symlinks require an ownership marker.
        if not (destination / ".omics-managed").is_file() and any(
            p.is_file() or p.is_symlink() for p in destination.rglob("*")
        ):
            raise ValueError(f"Refusing to replace unmanaged directory: {destination}")
        destination.replace(Path(scratch) / "previous")
    (stage / ".omics-managed").write_text("omics-workflows chromosome module\n", encoding="utf-8")
    stage.replace(destination)


@contextmanager
def assembly_scratch(destination, prefix=".joint-"):
    """在目标父目录创建暂存工作区并以字符串路径交给调用者。

    成功后清理临时文件；失败时移入failed_attempts保留现场并重新抛出原异常。
    """
    temporary = Path(tempfile.mkdtemp(prefix=prefix, dir=destination.parent))
    try:
        yield str(temporary)
    except BaseException:
        failed = destination.parent / "failed_attempts"
        failed.mkdir(exist_ok=True)
        retained = failed / temporary.name.lstrip(".")
        temporary.replace(retained)
        print(f"FAILED_ASSEMBLY_RETAINED {retained.resolve()}", flush=True)
        raise
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def pipeline(commands, target):
    """以参数列表连接多个外部进程，把末端二进制输出写入target。

    不经过shell插值；任何子进程失败均报错，包括SIGPIPE，并确保退出时回收所有子进程。
    """
    print("PIPELINE", json.dumps([[str(x) for x in c] for c in commands]), flush=True)
    processes = []
    try:
        with Path(target).open("wb") as out:
            previous = None
            for i, command in enumerate(commands):
                # Explicit process cleanup below owns the pipe/signal lifetime.
                # pylint: disable-next=consider-using-with
                p = subprocess.Popen(
                    [str(x) for x in command],
                    stdin=previous,
                    stdout=out if i == len(commands) - 1 else subprocess.PIPE,
                )
                if previous is not None:
                    previous.close()
                previous = p.stdout
                processes.append(p)
            codes = [p.wait() for p in processes]
            if any(codes):
                raise RuntimeError(f"Pipeline subprocess failure: {codes}")
    finally:
        for p in processes:
            if p.poll() is None:
                p.terminate()
        for p in processes:
            p.wait()
