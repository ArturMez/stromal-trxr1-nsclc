#!/usr/bin/env python
"""
A minimal random-access file object over HTTP byte ranges, for feeding to h5py.

Provides RangeFile, used by stage2_emtab13526_01_read_author_obs.py to read the obs of the two
97 GiB E-MTAB-13526 author h5ad files on the EBI FTP server without downloading them. Unlike fsspec's
HTTPFile, which can stall without a read timeout, it sets explicit connect and read timeouts, retries
each block, caches fetched blocks, and counts what it fetched so the cost of a "no download" read can
be reported. No random numbers are drawn.
"""
import os
import time
from collections import OrderedDict

import requests


class RangeFile:
    """A read-only, seekable file object whose reads are served by cached HTTP range requests."""
    def __init__(self, url, block_size=1 << 20, max_cache_blocks=4096,
                 timeout=(15, 60), retries=6, session=None, verbose=False):
        """Probe the remote size with a one-byte range request and set the block size, cache and timeouts."""
        self.url = url
        self.bs = block_size
        self.timeout = timeout
        self.retries = retries
        self.verbose = verbose
        self.s = session or requests.Session()
        self._cache = OrderedDict()
        self._max_blocks = max_cache_blocks
        self.pos = 0
        self.n_requests = 0
        self.bytes_fetched = 0
        self.wait_seconds = 0.0

        r = self._request_head()
        self.size = r

    def _request_head(self):
        """Return the remote file size from the Content-Range header of a one-byte request (with retries)."""
        last = None
        for a in range(self.retries):
            try:
                resp = self.s.get(self.url, headers={"Range": "bytes=0-0"},
                                  timeout=self.timeout, stream=True)
                resp.raise_for_status()
                cr = resp.headers.get("Content-Range", "")
                resp.close()
                if "/" not in cr:
                    raise RuntimeError(f"server did not honour Range: {resp.headers}")
                return int(cr.split("/")[-1])
            except Exception as e:  # noqa: BLE001
                last = e
                time.sleep(min(2 ** a, 20))
        raise RuntimeError(f"HEAD/range probe failed for {self.url}: {last}")

    def _fetch(self, bi):
        """Download block number bi (one range request, with retries) and return its bytes."""
        start = bi * self.bs
        end = min(start + self.bs, self.size) - 1
        last = None
        for a in range(self.retries):
            t0 = time.time()
            try:
                resp = self.s.get(self.url,
                                  headers={"Range": f"bytes={start}-{end}"},
                                  timeout=self.timeout)
                resp.raise_for_status()
                b = resp.content
                if len(b) != end - start + 1:
                    raise RuntimeError(f"short block {len(b)} != {end-start+1}")
                self.n_requests += 1
                self.bytes_fetched += len(b)
                self.wait_seconds += time.time() - t0
                return b
            except Exception as e:  # noqa: BLE001
                last = e
                self.wait_seconds += time.time() - t0
                if self.verbose:
                    print(f"      retry block {bi} ({a+1}/{self.retries}): {e}", flush=True)
                self.s = requests.Session()
                time.sleep(min(2 ** a, 20))
        raise RuntimeError(f"block {bi} failed after {self.retries} tries: {last}")

    def _block(self, bi):
        """Return block bi from the least-recently-used cache, fetching it when absent."""
        b = self._cache.get(bi)
        if b is None:
            b = self._fetch(bi)
            self._cache[bi] = b
            while len(self._cache) > self._max_blocks:
                self._cache.popitem(last=False)
        else:
            self._cache.move_to_end(bi)
        return b

    # ---- file protocol ----
    def readable(self):
        """File protocol: always True."""
        return True

    def seekable(self):
        """File protocol: always True."""
        return True

    def writable(self):
        """File protocol: always False."""
        return False

    def tell(self):
        """Return the current read position."""
        return self.pos

    def seek(self, off, whence=0):
        """Move the read position (whence 0, 1 or 2 as for a file)."""
        if whence == 0:
            self.pos = off
        elif whence == 1:
            self.pos += off
        elif whence == 2:
            self.pos = self.size + off
        return self.pos

    def read(self, n=-1):
        """Read n bytes (all remaining when n is negative) from the current position."""
        if n is None or n < 0:
            n = self.size - self.pos
        n = max(0, min(n, self.size - self.pos))
        if n == 0:
            return b""
        out = bytearray()
        p = self.pos
        remaining = n
        while remaining > 0:
            bi = p // self.bs
            off = p - bi * self.bs
            blk = self._block(bi)
            take = min(remaining, len(blk) - off)
            out += blk[off:off + take]
            p += take
            remaining -= take
        self.pos = p
        return bytes(out)

    def readinto(self, buf):
        """Fill buf from the current position and return the number of bytes read."""
        data = self.read(len(buf))
        buf[:len(data)] = data
        return len(data)

    def close(self):
        """Clear the cache and close the HTTP session."""
        self._cache.clear()
        try:
            self.s.close()
        except Exception:  # noqa: BLE001
            pass

    def __enter__(self):
        """Context manager entry; returns the file."""
        return self

    def __exit__(self, *a):
        """Context manager exit; closes the file."""
        self.close()

    def stats(self):
        """Return the number of requests, MiB fetched, network seconds and remote file size."""
        return {
            "requests": self.n_requests,
            "MiB_fetched": round(self.bytes_fetched / 1048576, 1),
            "network_seconds": round(self.wait_seconds, 1),
            "file_GiB": round(self.size / 1073741824, 2),
        }
