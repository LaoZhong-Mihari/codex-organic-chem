"""Literature search that returns something a chemist can judge: abstracts.

Titles alone cannot tell whether a paper is precedent for a transformation.
OpenAlex is queried first because it serves abstracts, citation counts and
open-access links without an API key; Crossref fills in when OpenAlex is
unreachable or returns zero hits. Only the query text is sent.

An abstract is still not a procedure. Scope, yields and failures live in the
paper and its supporting information; cite a DOI as precedent only for what the
abstract (or the paper, once read) actually states.
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import OrderedDict
from typing import Any


CROSSREF_WORKS_URL = "https://api.crossref.org/works"
OPENALEX_WORKS_URL = "https://api.openalex.org/works"
USER_AGENT = "codex-organic-chem/1.1 (local research assistant)"
ABSTRACT_LIMIT = 1500
MAX_CACHE_ENTRIES = 128

_RESPONSE_CACHE: OrderedDict[str, dict[str, Any]] = OrderedDict()


def _first(value: Any) -> Any:
    if isinstance(value, list):
        return value[0] if value else None
    return value


def _user_agent() -> str:
    mailto = os.environ.get("CODEX_CHEM_MAILTO", "").strip()
    if mailto:
        return f"{USER_AGENT} (mailto:{mailto})"
    return USER_AGENT


def _get_json(url: str, timeout_s: int) -> dict[str, Any]:
    if url in _RESPONSE_CACHE:
        _RESPONSE_CACHE.move_to_end(url)
        return _RESPONSE_CACHE[url]
    headers = {"Accept": "application/json", "User-Agent": _user_agent()}
    delays = (0.5, 1.2)
    last_exc: Exception | None = None
    for attempt in range(len(delays) + 1):
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=timeout_s) as response:
                data = json.loads(response.read().decode("utf-8"))
            _RESPONSE_CACHE[url] = data
            while len(_RESPONSE_CACHE) > MAX_CACHE_ENTRIES:
                _RESPONSE_CACHE.popitem(last=False)
            return data
        except urllib.error.HTTPError as exc:
            last_exc = exc
            if exc.code in {429, 502, 503, 504} and attempt < len(delays):
                time.sleep(delays[attempt])
                continue
            raise
        except Exception as exc:
            last_exc = exc
            raise
    assert last_exc is not None
    raise last_exc


def _normalize_doi(doi: str | None) -> str | None:
    if not doi:
        return None
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", doi.strip(), flags=re.IGNORECASE).lower()


def _abstract_from_inverted_index(index: dict[str, list[int]] | None) -> str | None:
    if not index:
        return None
    words = sorted((position, word) for word, positions in index.items() for position in positions)
    text = " ".join(word for _, word in words)
    return text if len(text) <= ABSTRACT_LIMIT else text[:ABSTRACT_LIMIT].rsplit(" ", 1)[0] + " ..."


def _strip_jats(text: str | None) -> str | None:
    if not text:
        return None
    clean = re.sub(r"<[^>]+>", " ", text)
    clean = re.sub(r"\s+", " ", clean).strip()
    return clean if len(clean) <= ABSTRACT_LIMIT else clean[:ABSTRACT_LIMIT].rsplit(" ", 1)[0] + " ..."


def search_openalex(
    query: str, rows: int = 5, timeout_s: int = 15, year_from: int | None = None
) -> dict[str, Any]:
    try:
        per_page = max(1, min(int(rows), 25))
        yr = int(year_from) if year_from is not None else None
    except (TypeError, ValueError) as exc:
        return {
            "source": "OpenAlex",
            "query": query,
            "status": "error",
            "url": OPENALEX_WORKS_URL,
            "warnings": [f"Invalid numeric parameter: {exc}"],
            "items": [],
        }
    params: dict[str, Any] = {
        "search": query,
        "per-page": per_page,
        "select": "id,doi,title,publication_year,primary_location,abstract_inverted_index,cited_by_count,type,open_access",
    }
    if yr:
        params["filter"] = f"from_publication_date:{yr}-01-01"
    mailto = os.environ.get("CODEX_CHEM_MAILTO", "").strip()
    if mailto:
        params["mailto"] = mailto
    api_key = os.environ.get("OPENALEX_API_KEY", "").strip()
    if api_key:
        params["api_key"] = api_key
    url = f"{OPENALEX_WORKS_URL}?{urllib.parse.urlencode(params)}"
    try:
        payload = _get_json(url, timeout_s)
    except Exception as exc:
        return {
            "source": "OpenAlex",
            "query": query,
            "status": "error",
            "url": url,
            "warnings": [f"OpenAlex request failed: {exc}"],
            "items": [],
        }
    items = []
    for work in payload.get("results", []):
        doi = _normalize_doi(work.get("doi"))
        location = work.get("primary_location") or {}
        source = location.get("source") or {}
        items.append(
            {
                "title": work.get("title"),
                "journal": source.get("display_name"),
                "year": work.get("publication_year"),
                "doi": doi,
                "url": f"https://doi.org/{doi}" if doi else work.get("id"),
                "type": work.get("type"),
                "cited_by": work.get("cited_by_count"),
                "open_access_url": (work.get("open_access") or {}).get("oa_url"),
                "abstract": _abstract_from_inverted_index(work.get("abstract_inverted_index")),
                "source": "OpenAlex",
            }
        )
    return {"source": "OpenAlex", "query": query, "status": "ok", "url": url, "items": items, "warnings": []}


def search_crossref(
    query: str, rows: int = 5, timeout_s: int = 12, year_from: int | None = None
) -> dict[str, Any]:
    try:
        limit = max(1, min(int(rows), 10))
        yr = int(year_from) if year_from is not None else None
    except (TypeError, ValueError) as exc:
        return {
            "source": "Crossref",
            "query": query,
            "status": "error",
            "url": CROSSREF_WORKS_URL,
            "warnings": [f"Invalid numeric parameter: {exc}"],
            "items": [],
        }
    query_params: dict[str, Any] = {
        "query.bibliographic": query,
        "rows": limit,
        "select": "DOI,title,container-title,published-print,published-online,issued,URL,score,type,abstract",
    }
    if yr:
        query_params["filter"] = f"from-pub-date:{yr}"
    mailto = os.environ.get("CODEX_CHEM_MAILTO", "").strip()
    if mailto:
        query_params["mailto"] = mailto
    url = f"{CROSSREF_WORKS_URL}?{urllib.parse.urlencode(query_params)}"
    try:
        payload = _get_json(url, timeout_s)
    except Exception as exc:
        return {
            "source": "Crossref",
            "query": query,
            "status": "error",
            "url": url,
            "warnings": [f"Crossref request failed: {exc}"],
            "items": [],
        }
    items = []
    for item in payload.get("message", {}).get("items", []):
        issued = item.get("published-print") or item.get("published-online") or item.get("issued") or {}
        date_parts = issued.get("date-parts", [[]])
        year = date_parts[0][0] if date_parts and date_parts[0] else None
        if yr and (year is None or year < yr):
            continue
        doi = _normalize_doi(item.get("DOI"))
        items.append(
            {
                "title": _first(item.get("title")),
                "journal": _first(item.get("container-title")),
                "year": year,
                "doi": doi,
                "url": item.get("URL") or (f"https://doi.org/{doi}" if doi else None),
                "type": item.get("type"),
                "score": item.get("score"),
                "abstract": _strip_jats(item.get("abstract")),
                "source": "Crossref",
            }
        )
    return {
        "source": "Crossref",
        "query": query,
        "status": "ok",
        "url": url,
        "items": items,
        "warnings": [],
    }


def literature_search(
    query: str, rows: int = 5, source: str = "auto", year_from: int | None = None
) -> dict[str, Any]:
    """Search OpenAlex (abstracts) with Crossref as fallback when OpenAlex fails or finds nothing."""
    if source not in {"auto", "openalex", "crossref"}:
        return {"status": "error", "query": query, "items": [], "warnings": [f"Unknown source {source!r}."]}
    try:
        target_rows = max(1, int(rows))
        yr = int(year_from) if year_from is not None else None
    except (TypeError, ValueError) as exc:
        return {"status": "error", "query": query, "items": [], "warnings": [f"Invalid parameter: {exc}"]}
    warnings: list[str] = []
    items: list[dict[str, Any]] = []
    sources: list[str] = []
    urls: list[str] = []
    if source in {"auto", "openalex"}:
        primary = search_openalex(query, rows=target_rows, year_from=yr)
        warnings.extend(primary["warnings"])
        urls.append(primary["url"])
        if primary["status"] == "ok":
            sources.append("OpenAlex")
            items.extend(primary["items"])
    if source == "crossref" or (source == "auto" and len(items) < target_rows):
        secondary = search_crossref(query, rows=target_rows, year_from=yr)
        warnings.extend(secondary["warnings"])
        urls.append(secondary["url"])
        if secondary["status"] == "ok":
            sources.append("Crossref")
            seen = {item["doi"] for item in items if item.get("doi")}
            had_primary = bool(items)
            for item in secondary["items"]:
                if item.get("doi") and item["doi"] in seen:
                    continue
                if yr and (item.get("year") is None or item["year"] < yr):
                    continue
                if had_primary and item.get("score") is not None and float(item["score"]) < 15.0:
                    continue
                items.append(item)
    items = items[:target_rows]
    without_abstract = sum(1 for item in items if not item.get("abstract"))
    if items and without_abstract:
        warnings.append(
            f"{without_abstract} of {len(items)} records have no abstract; a title alone does not establish precedent."
        )
    return {
        "source": "+".join(sources) or "none",
        "query": query,
        "status": "ok" if sources else "error",
        "url": urls[0] if urls else None,
        "items": items,
        "warnings": warnings,
        "evidence_note": (
            "Abstract-level evidence. Grade each hit as exact precedent, close analogue, same reaction class "
            "or review, and read the paper/SI before relying on scope, selectivity or conditions."
        ),
    }
