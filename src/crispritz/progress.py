""" """

from .verbosity import VERBOSITY_LVL

from typing import Iterable
from threading import Lock
from contextlib import contextmanager
from tqdm import tqdm


class ThreadSafeProgressBar:

    def __init__(self, pbar: tqdm) -> None:
        self._pbar = pbar
        self._lock = Lock()

    def update(self, n: int = 1) -> None:
        with self._lock:
            self._pbar.update(n)

    def set_description(self, desc: str) -> None:
        with self._lock:
            self._pbar.set_description(desc)

    def close(self) -> None:
        with self._lock:
            self._pbar.close()


class DummyProgressBar:

    def update(self, n: int = 1) -> None:
        pass

    def set_description(self, desc: str) -> None:
        pass

    def close(self) -> None:
        pass


def progress_bar(iterable: Iterable, desc: str, verbosity: int) -> Iterable:
    if verbosity in VERBOSITY_LVL[1:3]:
        return tqdm(
            iterable,
            desc=desc,
            ascii="░█",
            bar_format="{desc}: |{bar:100}| {percentage:3.0f}%",
        )
    return iterable


@contextmanager
def progress_bar_parallel(total: int, desc: str, verbosity: int):
    if verbosity in VERBOSITY_LVL[1:3]:
        pbar = tqdm(
            total=total,
            desc=desc,
            ascii="░█",
            initial=0,
            bar_format="{desc}: |{bar:100}| {percentage:3.0f}%",
        )
        try:
            yield ThreadSafeProgressBar(pbar)
        finally:
            pbar.close()
    else:
        yield DummyProgressBar()
