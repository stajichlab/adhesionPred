import os
from concurrent.futures import ProcessPoolExecutor

import pytest

from cellsurface_sorting_hat.cache import (
    CacheError,
    ModuleCache,
    identity_key,
    read_verified,
    write_atomic,
)
from cellsurface_sorting_hat.status import ModuleIdentity

IDENT = ModuleIdentity("pfam", "1", "p", "a")


def test_identity_key_changes_with_every_field():
    base = identity_key(IDENT, {"hmmer": "3.4"})
    assert base == identity_key(IDENT, {"hmmer": "3.4"})
    for changed in (
        ModuleIdentity("pfam", "2", "p", "a"),
        ModuleIdentity("pfam", "1", "q", "a"),
        ModuleIdentity("pfam", "1", "p", "b"),
        ModuleIdentity("other", "1", "p", "a"),
    ):
        assert identity_key(changed, {"hmmer": "3.4"}) != base
    assert identity_key(IDENT, {"hmmer": "3.5"}) != base


def test_write_atomic_leaves_data_and_checksum_and_no_temporary_file(tmp_path):
    path = tmp_path / "d" / "f.bin"
    write_atomic(path, b"abc")
    assert read_verified(path) == b"abc"
    assert sorted(p.name for p in path.parent.iterdir()) == ["f.bin", "f.bin.sha256"]


def test_read_verified_refuses_a_missing_or_wrong_checksum(tmp_path):
    path = tmp_path / "f.bin"
    path.write_bytes(b"abc")
    with pytest.raises(CacheError, match="no checksum"):
        read_verified(path)
    write_atomic(path, b"abc")
    path.write_bytes(b"abd")
    with pytest.raises(CacheError, match="does not match"):
        read_verified(path)


def test_module_cache_merges_rows_by_sha256(tmp_path):
    cache = ModuleCache(tmp_path, "k1")
    assert cache.load() == {}
    cache.update([{"sha256": "aa", "call": "called"}, {"sha256": "bb", "call": "not_called"}])
    cache.update([{"sha256": "bb", "call": "called"}, {"sha256": "cc", "call": "called", "x": "1"}])
    table = ModuleCache(tmp_path, "k1").load()
    assert set(table) == {"aa", "bb", "cc"}
    assert table["bb"]["call"] == "called"
    assert table["cc"]["x"] == "1"
    assert ModuleCache(tmp_path, "k2").load() == {}  # another identity key: another table


def test_a_damaged_table_is_a_cache_miss_and_the_next_update_rebuilds_it(tmp_path):
    cache = ModuleCache(tmp_path, "k1")
    cache.update([{"sha256": "aa", "call": "called"}])
    cache.path.write_bytes(b"not a gzip file")  # the sidecar no longer matches
    assert cache.load() == {}
    cache.update([{"sha256": "bb", "call": "called"}])
    assert set(cache.load()) == {"bb"}


def _worker(args):
    directory, worker = args
    cache = ModuleCache(directory, "shared")
    for i in range(15):
        cache.update([{"sha256": f"w{worker}-{i}", "call": "called"}])
        cache.load()  # readers run while other processes write
    return worker


def test_parallel_updates_lose_no_rows_and_leave_a_readable_table(tmp_path):
    with ProcessPoolExecutor(max_workers=4) as pool:
        list(pool.map(_worker, [(str(tmp_path), w) for w in range(4)]))
    table = ModuleCache(tmp_path, "shared").load()
    assert len(table) == 60
    assert read_verified(ModuleCache(tmp_path, "shared").path)  # data and checksum agree


def test_a_reader_can_read_when_the_lock_file_cannot_be_opened(tmp_path):
    cache = ModuleCache(tmp_path, "k1")
    cache.update([{"sha256": "aa", "call": "called"}])
    lock = cache._lock
    lock.chmod(0o444)
    try:
        if os.access(lock, os.W_OK):
            pytest.skip("this user can write to a 0444 file (for example root)")
        assert set(cache.load()) == {"aa"}
    finally:
        lock.chmod(0o644)
