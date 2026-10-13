"""Running one computation in a child process with memory and wall-clock limits.

Used by the long table scripts so that a single hard field (an LLL or MILP that
blows up, or a FLINT abort on a failed allocation) cannot stall or kill a run.
"""

from __future__ import annotations

import multiprocessing as mp
import queue
import time


def _child(q, target, args, mem_bytes):
    import resource

    if mem_bytes:
        resource.setrlimit(resource.RLIMIT_AS, (mem_bytes, mem_bytes))
    try:
        q.put(("ok", target(*args)))
    except MemoryError:
        q.put(("error", "memory limit"))
    except Exception as e:  # report, do not lose the run
        q.put(("error", f"{type(e).__name__}: {e}"))


def run_limited(target, args=(), mem_bytes: int = 5 << 30, wall: float = 1800.0):
    """``(target(*args), None)``, or ``(None, reason)`` if the child failed or ran out of time.

    The result must be picklable.  A child killed by a signal (e.g. ``-6``: FLINT
    aborts when an allocation fails under the memory limit) is reported by its
    exit code.
    """
    q = mp.Queue()
    p = mp.Process(target=_child, args=(q, target, args, mem_bytes))
    p.start()
    deadline = time.time() + wall
    while True:
        try:
            status, val = q.get(timeout=1.0)
        except queue.Empty:
            if not p.is_alive():
                try:
                    status, val = q.get(timeout=1.0)
                except queue.Empty:
                    p.join()
                    return None, f"child exit code {p.exitcode}"
            elif time.time() > deadline:
                p.kill()
                p.join()
                return None, "wall-clock limit"
            else:
                continue
        p.join()
        return (val, None) if status == "ok" else (None, val)
