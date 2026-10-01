import gzip
import hashlib
import http.client
import io
import urllib.error

import manifest
import paths
import pytest
from conftest import load_script

GOOD = gzip.compress(b"!gaf-version: 2.2\nSGD\tS1\n", mtime=0)
GOOD_SHA = hashlib.sha256(GOOD).hexdigest()


class FakeResponse(io.BytesIO):
    def __init__(self, body, headers):
        super().__init__(body)
        self.headers = headers

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def fake_opener(pages):
    """pages: url -> (body, headers). Records every requested URL."""
    calls = []

    def opener(request, timeout=None):
        calls.append(request.full_url)
        body, headers = pages[request.full_url]
        return FakeResponse(body, headers)

    opener.calls = calls
    return opener


def write_manifest(path, rows):
    full = [{c: "" for c in manifest.MANIFEST_COLUMNS} | r for r in rows]
    manifest.write_manifest(path, full)
    return path


def test_check_payload_refuses_html_truncated_and_empty(tmp_path):
    html = tmp_path / "a.gaf.gz"
    html.write_bytes(b"<!DOCTYPE html><html><body>403 Forbidden</body></html>")
    with pytest.raises(manifest.DownloadError, match="HTML"):
        manifest.check_payload(html)
    html_gz = tmp_path / "b.gaf.gz"
    html_gz.write_bytes(gzip.compress(b"<html>error</html>"))
    with pytest.raises(manifest.DownloadError, match="HTML"):
        manifest.check_payload(html_gz)
    cut = tmp_path / "c.gaf.gz"
    noise = b"".join(hashlib.sha256(str(i).encode()).digest() for i in range(3000))
    cut.write_bytes(gzip.compress(noise)[:50000])
    with pytest.raises(manifest.DownloadError, match="truncated"):
        manifest.check_payload(cut)
    empty = tmp_path / "d.gaf.gz"
    empty.write_bytes(b"")
    with pytest.raises(manifest.DownloadError, match="empty"):
        manifest.check_payload(empty)


def test_strict_hash_mismatch_stops_and_keeps_no_file(tmp_path):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/sgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "sgd.gaf.gz", "kind": "http", "url": url, "sha256": "0" * 64, "mode": "strict"}],
    )
    opener = fake_opener({url: (GOOD, {"Last-Modified": "Thu, 28 May 2026 21:18:22 GMT"})})
    dest = tmp_path / "dl"
    assert fetch.main(["--manifest", str(mpath), "--dest", str(dest)], opener=opener) == 2
    assert not (dest / "sgd.gaf.gz").exists()
    assert not (dest / "sgd.gaf.gz.part").exists()


def test_update_manifest_accepts_new_hash(tmp_path):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/sgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "sgd.gaf.gz", "kind": "http", "url": url, "sha256": "0" * 64, "mode": "strict"}],
    )
    opener = fake_opener({url: (GOOD, {"Last-Modified": "Thu, 28 May 2026 21:18:22 GMT"})})
    argv = ["--manifest", str(mpath), "--dest", str(tmp_path / "dl"), "--update-manifest"]
    assert fetch.main(argv, opener=opener) == 0
    row = manifest.read_manifest(mpath)[0]
    assert row["sha256"] == GOOD_SHA
    assert row["last_modified"] == "Thu, 28 May 2026 21:18:22 GMT"
    assert (tmp_path / "dl" / "fetch_log.tsv").read_text().count("sgd.gaf.gz") == 2


def test_html_error_page_stops_even_with_update(tmp_path):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/aspgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "aspgd.gaf.gz", "kind": "http", "url": url, "sha256": "", "mode": "record"}],
    )
    opener = fake_opener({url: (b"<html><body>403</body></html>", {})})
    argv = ["--manifest", str(mpath), "--dest", str(tmp_path / "dl"), "--update-manifest"]
    assert fetch.main(argv, opener=opener) == 2
    assert not (tmp_path / "dl" / "aspgd.gaf.gz").exists()


def test_cached_strict_file_is_not_downloaded_again(tmp_path):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/sgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "sgd.gaf.gz", "kind": "http", "url": url, "sha256": GOOD_SHA, "mode": "strict"}],
    )
    (tmp_path / "dl").mkdir()
    (tmp_path / "dl" / "sgd.gaf.gz").write_bytes(GOOD)
    opener = fake_opener({})
    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 0
    )
    assert opener.calls == []


def test_uniprot_pages_follow_link_header(tmp_path):
    fetch = load_script("00_fetch_inputs")
    first = "https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP1&format=fasta&size=500"
    second = "https://rest.uniprot.org/uniprotkb/search?cursor=abc"
    pages = {
        first: (
            b">sp|P1|A_B\nMK\n",
            {"Link": f'<{second}>; rel="next"', "X-UniProt-Release": "2026_03"},
        ),
        second: (b">tr|Q2|C_D\nMA\n", {}),
    }
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "UP1.fasta.gz", "kind": "uniprot_fasta", "url": first, "mode": "record"}],
    )
    opener = fake_opener(pages)
    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 0
    )
    assert opener.calls == [first, second]
    assert gzip.decompress((tmp_path / "dl" / "UP1.fasta.gz").read_bytes()).count(b">") == 2


def test_http_403_stops_with_exit_2(tmp_path, capsys):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/aspgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "aspgd.gaf.gz", "kind": "http", "url": url, "sha256": "", "mode": "strict"}],
    )

    def opener(request, timeout=None):
        raise urllib.error.HTTPError(request.full_url, 403, "Forbidden", {}, None)

    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 2
    )
    assert "STOP:" in capsys.readouterr().err
    assert list((tmp_path / "dl").iterdir()) == []


def test_url_error_mid_download_leaves_no_part_file(tmp_path):
    fetch = load_script("00_fetch_inputs")
    first = "https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP1&format=fasta&size=500"
    second = "https://rest.uniprot.org/uniprotkb/search?cursor=abc"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "UP1.fasta.gz", "kind": "uniprot_fasta", "url": first, "mode": "record"}],
    )

    def opener(request, timeout=None):
        if request.full_url == first:
            return FakeResponse(b">sp|P1|A_B\nMK\n", {"Link": f'<{second}>; rel="next"'})
        raise urllib.error.URLError("connection reset")

    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 2
    )
    assert list((tmp_path / "dl").iterdir()) == []


def test_committed_manifest_covers_species_table():
    rows = {r["file"]: r for r in manifest.read_manifest(paths.STEP1_DIR / "manifest.tsv")}
    import csv

    with open(paths.STEP1_DIR / "species.tsv", encoding="utf-8") as handle:
        species = list(csv.DictReader(handle, delimiter="\t"))
    assert "go-basic.obo" in rows
    for sp in species:
        assert rows[sp["gaf_file"]]["mode"] == "strict"
        assert len(rows[sp["gaf_file"]]["sha256"]) == 64
        assert sp["fasta_file"] in rows


class FailingRead(FakeResponse):
    def __init__(self, exc):
        super().__init__(b"", {})
        self.exc = exc

    def read(self, *args):
        raise self.exc


def _single_row(tmp_path, kind="http", mode="strict", sha=""):
    url = "https://example.org/x.gaf.gz"
    if kind == "uniprot_fasta":
        url = "https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP1&format=fasta&size=500"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "x.gaf.gz", "kind": kind, "url": url, "sha256": sha, "mode": mode}],
    )
    return mpath, url


@pytest.mark.parametrize(
    "exc",
    [ConnectionResetError("reset by peer"), http.client.IncompleteRead(b"ab", 10), TimeoutError()],
)
def test_read_error_gives_stop_and_no_part_file(tmp_path, capsys, exc):
    fetch = load_script("00_fetch_inputs")
    mpath, url = _single_row(tmp_path)

    def opener(request, timeout=None):
        return FailingRead(exc)

    dest = tmp_path / "dl"
    assert fetch.main(["--manifest", str(mpath), "--dest", str(dest)], opener=opener) == 2
    err = capsys.readouterr().err
    assert "STOP:" in err
    assert url in err
    assert list(dest.iterdir()) == []


def test_error_on_second_uniprot_page_leaves_no_output(tmp_path, capsys):
    fetch = load_script("00_fetch_inputs")
    mpath, first = _single_row(tmp_path, kind="uniprot_fasta", mode="record")
    second = "https://rest.uniprot.org/uniprotkb/search?cursor=abc"

    def opener(request, timeout=None):
        if request.full_url == first:
            return FakeResponse(b">sp|P1|A_B\nMK\n", {"Link": f'<{second}>; rel="next"'})
        return FailingRead(ConnectionResetError("reset"))

    dest = tmp_path / "dl"
    assert fetch.main(["--manifest", str(mpath), "--dest", str(dest)], opener=opener) == 2
    err = capsys.readouterr().err
    assert "STOP:" in err
    assert first in err
    assert list(dest.iterdir()) == []


def test_verify_strict_mismatch_names_file_and_hashes(tmp_path):
    path = tmp_path / "f.gaf.gz"
    path.write_bytes(GOOD)
    row = {"mode": "strict", "sha256": "ab" * 32}
    with pytest.raises(manifest.DownloadError) as info:
        manifest.verify_against_manifest(path, row)
    text = str(info.value)
    assert "f.gaf.gz" in text and GOOD_SHA[:12] in text and "abababababab" in text


def test_verify_strict_match_returns_sha(tmp_path):
    path = tmp_path / "f.gaf.gz"
    path.write_bytes(GOOD)
    assert (
        manifest.verify_against_manifest(path, {"mode": "strict", "sha256": GOOD_SHA}) == GOOD_SHA
    )


def test_verify_record_mode_returns_new_hash(tmp_path):
    path = tmp_path / "f.gaf.gz"
    path.write_bytes(GOOD)
    row = {"mode": "record", "sha256": "0" * 64}
    assert manifest.verify_against_manifest(path, row) == GOOD_SHA


def test_verify_uppercase_manifest_hash_matches(tmp_path):
    path = tmp_path / "f.gaf.gz"
    path.write_bytes(GOOD)
    row = {"mode": "strict", "sha256": GOOD_SHA.upper()}
    assert manifest.verify_against_manifest(path, row) == GOOD_SHA


def test_uppercase_hash_in_main_is_cached(tmp_path):
    fetch = load_script("00_fetch_inputs")
    mpath, _ = _single_row(tmp_path, sha=GOOD_SHA.upper())
    (tmp_path / "dl").mkdir()
    (tmp_path / "dl" / "x.gaf.gz").write_bytes(GOOD)
    opener = fake_opener({})
    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 0
    )
    assert opener.calls == []


def test_failure_message_names_url_and_file(tmp_path, capsys):
    fetch = load_script("00_fetch_inputs")
    mpath, url = _single_row(tmp_path, sha="0" * 64)
    opener = fake_opener({url: (GOOD, {})})
    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 2
    )
    err = capsys.readouterr().err
    assert f"FAILED {url} -> x.gaf.gz" in err


def test_only_with_unknown_file_stops(tmp_path, capsys):
    fetch = load_script("00_fetch_inputs")
    mpath, _ = _single_row(tmp_path)
    argv = ["--manifest", str(mpath), "--dest", str(tmp_path / "dl"), "--only", "nope.gz"]
    assert fetch.main(argv, opener=fake_opener({})) == 2
    assert "STOP: unknown file nope.gz" in capsys.readouterr().err


def test_write_manifest_failure_keeps_old_file(tmp_path):
    path = write_manifest(tmp_path / "m.tsv", [{"file": "a", "mode": "strict"}])
    before = path.read_bytes()
    bad = [{c: "" for c in manifest.MANIFEST_COLUMNS} | {"file": "b", "bogus": "x"}]
    with pytest.raises(ValueError):
        manifest.write_manifest(path, bad)
    assert path.read_bytes() == before
    assert [p.name for p in tmp_path.iterdir()] == ["m.tsv"]
