"""Regression tests for Phase 5 literature search, packaging, and Python >=3.12 support."""

from __future__ import annotations

import io
import json
import urllib.error
from unittest.mock import patch

from codex_organic_chem import external, literature


def test_literature_search_polite_mailto_caching_retry_and_auto_fallback(monkeypatch):
    literature._RESPONSE_CACHE.clear()
    monkeypatch.setenv("CODEX_CHEM_MAILTO", "chemist@example.org")
    monkeypatch.setenv("OPENALEX_API_KEY", "oa-test-key")

    seen_requests = []
    attempts = {"count": 0}

    def fake_urlopen(req, timeout=15):
        seen_requests.append((req.full_url, dict(req.header_items())))
        attempts["count"] += 1
        if attempts["count"] == 1:
            raise urllib.error.HTTPError(req.full_url, 429, "Too Many Requests", {}, None)
        payload = {
            "results": [
                {
                    "id": "https://openalex.org/W1",
                    "doi": "https://doi.org/10.1021/jacs.0c00001",
                    "title": "Sharpless Asymmetric Epoxidation",
                    "publication_year": 2021,
                    "primary_location": {"source": {"display_name": "J. Am. Chem. Soc."}},
                    "abstract_inverted_index": {"Concerted": [0], "syn": [1], "epoxidation.": [2]},
                    "cited_by_count": 42,
                    "type": "article",
                    "open_access": {"oa_url": None},
                }
            ]
        }
        return io.BytesIO(json.dumps(payload).encode("utf-8"))

    with (
        patch("urllib.request.urlopen", side_effect=fake_urlopen),
        patch("time.sleep"),
    ):
        res = literature.literature_search("Sharpless epoxidation", rows=1, source="auto", year_from=2020)
        # Second identical call hits the in-memory LRU cache without new HTTP requests
        res_cached = literature.literature_search("Sharpless epoxidation", rows=1, source="auto", year_from=2020)

    assert res["status"] == "ok"
    assert res["source"] == "OpenAlex"
    assert len(res["items"]) == 1
    assert res["items"][0]["abstract"] == "Concerted syn epoxidation."
    assert res_cached["items"] == res["items"]
    assert attempts["count"] == 2  # 1 failed (429) + 1 succeeded; 3rd call served from cache
    url0, headers0 = seen_requests[0]
    assert "mailto=chemist%40example.org" in url0
    assert "api_key=oa-test-key" in url0
    assert "mailto:chemist@example.org" in headers0.get("User-agent", "")


def test_crossref_applies_year_filter_and_handles_bad_parameters():
    literature._RESPONSE_CACHE.clear()
    assert literature.literature_search("aldol", rows="not-an-int")["status"] == "error"

    def fake_crossref(req, timeout=12):
        assert "filter=from-pub-date%3A2022" in req.full_url
        payload = {
            "message": {
                "items": [
                    {"DOI": "10.1000/old", "title": ["Old"], "issued": {"date-parts": [[2019]]}},
                    {"DOI": "10.1000/ undated", "title": ["Undated"], "issued": {"date-parts": [[]]}},
                    {"DOI": "10.1000/new", "title": ["New"], "issued": {"date-parts": [[2023]]}, "abstract": "<jats:p>Clean.</jats:p>"},
                ]
            }
        }
        return io.BytesIO(json.dumps(payload).encode("utf-8"))

    with patch("urllib.request.urlopen", side_effect=fake_crossref):
        res = literature.literature_search("aldol", rows=5, source="crossref", year_from=2022)
    assert res["status"] == "ok"
    assert [it["year"] for it in res["items"]] == [2023]
    assert res["items"][0]["abstract"] == "Clean."


def test_doctor_accepts_python_312_and_313(monkeypatch):
    monkeypatch.setattr("sys.version_info", (3, 13, 2, "final", 0))
    env = external.python_environment()
    assert env["python"]["status"] == "available"
