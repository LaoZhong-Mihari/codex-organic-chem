from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import Any


CROSSREF_WORKS_URL = "https://api.crossref.org/works"


def _first(value: Any) -> Any:
    if isinstance(value, list):
        return value[0] if value else None
    return value


def search_crossref(query: str, rows: int = 5, timeout_s: int = 12) -> dict[str, Any]:
    params = urllib.parse.urlencode(
        {
            "query.bibliographic": query,
            "rows": max(1, min(int(rows), 10)),
            "select": "DOI,title,container-title,published-print,published-online,issued,URL,score,type",
        }
    )
    url = f"{CROSSREF_WORKS_URL}?{params}"
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "codex-organic-chem/0.1.0 (local research assistant)",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout_s) as response:
            payload = json.loads(response.read().decode("utf-8"))
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
        doi = item.get("DOI")
        items.append(
            {
                "title": _first(item.get("title")),
                "journal": _first(item.get("container-title")),
                "year": year,
                "doi": doi,
                "url": item.get("URL") or (f"https://doi.org/{doi}" if doi else None),
                "type": item.get("type"),
                "score": item.get("score"),
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


def literature_search(query: str, rows: int = 5) -> dict[str, Any]:
    return search_crossref(query=query, rows=rows)

