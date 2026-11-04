"""Keep failed native mapping/scaffolding work without changing assembly helpers."""
from contextlib import contextmanager
from pathlib import Path
import shutil
import tempfile


@contextmanager
def native_scratch(destination, prefix):
    temporary = Path(tempfile.mkdtemp(prefix=prefix, dir=destination.parent))
    try:
        yield str(temporary)
    except BaseException:
        failed = destination.parent/'failed_attempts'
        failed.mkdir(exist_ok=True)
        retained = failed/temporary.name.lstrip('.')
        temporary.replace(retained)
        print(f'FAILED_NATIVE_WORK_RETAINED {retained.resolve()}', flush=True)
        raise
    finally:
        if temporary.exists(): shutil.rmtree(temporary)
