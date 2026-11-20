"""A guided try-out of the review page with planted mistakes.

Opens two sessions: three structures "read from a scheme", one with a wedge
misread against its source image, and an eight-step route in which one
intermediate is drawn with the wrong configuration, so the audit objects until
the reviewer clicks the centre. Nothing here touches the user's own work.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .structure_review import home_dir, open_review

TRUE_PROLINE = "COC(=O)[C@@H]1C[C@H](O)CN1C(=O)OC(C)(C)C"  # cis: what the "source" shows
MISREAD_PROLINE = "COC(=O)[C@@H]1C[C@@H](O)CN1C(=O)OC(C)(C)C"  # trans: the planted misreading

_TEXT = {
    "en": {
        "s_title": "Demo: structures read from a scheme",
        "s_intro": (
            "Practice session with a planted mistake. Compare compound 2 with its source image: the C4 hydroxyl "
            "was read with the wrong wedge. Click that stereocentre in the drawing to flip it, then confirm. "
            "Compound 3 has an open stereocentre and an open double bond; set them with the R/S and E/Z buttons."
        ),
        "l1": "1 · Roche ester", "l2": "2 · N-Boc hydroxyproline methyl ester", "l3": "3 · allylic alcohol",
        "src": "Scheme 1",
        "n2": "C4-OH read as a hashed bond (trans to the ester). The bond is faint in the scan.",
        "q2": "Is the C4 hydroxyl a wedge or a hash in your drawing?",
        "n3": "No wedge and no alkene geometry is visible in the source; left undefined.",
        "alt": "C4 epimer (cis)",
        "r_title": "Demo: hydroxyproline to a 4-acetamido acid (planted error in compound 7)",
        "r_intro": (
            "Practice session. Step 5 is an SN2 azide displacement, but compound 7 is drawn with retention, so the "
            "audit objects at steps 5 and 6. Click compound 7, click its C4 centre to flip it, apply, and watch "
            "the audit clear. Then mark steps, answer the audit prompts and approve or send back."
        ),
    },
    "zh": {
        "s_title": "演示：从路线图读出的结构",
        "s_intro": (
            "这是一个带有预设错误的练习。请把化合物 2 与来源图对照：C4 羟基的楔形键被读反了。"
            "在结构图上点击这个立体中心即可翻转，然后确认。化合物 3 有一个未定义的立体中心和一个未定义的双键，"
            "可以用 R/S、E/Z 按钮指定。"
        ),
        "l1": "1 · Roche 酯", "l2": "2 · N-Boc 羟脯氨酸甲酯", "l3": "3 · 烯丙醇",
        "src": "Scheme 1",
        "n2": "C4-OH 读作虚楔键（与酯基反式）。扫描件里这根键比较淡。",
        "q2": "你的原图里 C4 羟基是实楔键还是虚楔键？",
        "n3": "来源图中看不到楔形键和双键构型，暂未定义。",
        "alt": "C4 差向异构体（顺式）",
        "r_title": "演示：羟脯氨酸到 4-乙酰氨基酸（化合物 7 有预设错误）",
        "r_intro": (
            "练习用。第 5 步是叠氮的 SN2 取代，但化合物 7 被画成了构型保持，所以核查在第 5 步和第 6 步报问题。"
            "点击化合物 7，再点它的 C4 中心翻转，应用后可以看到核查结果恢复正常。然后给各步打勾、回答核查提示，"
            "最后通过或退回。"
        ),
    },
}


def _source_image(directory: Path) -> str | None:
    """Draw the true structure the way a scheme would show it."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D

        mol = Chem.MolFromSmiles(TRUE_PROLINE)
        drawer = rdMolDraw2D.MolDraw2DCairo(560, 400)
        drawer.drawOptions().bondLineWidth = 2
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "scheme1_compound2.png"
        path.write_bytes(drawer.GetDrawingText())
        return str(path)
    except Exception:
        return None


def demo_route() -> dict[str, Any]:
    route = json.loads((Path(__file__).resolve().parent / "review_ui" / "demo_route.json").read_text(encoding="utf-8"))
    # planted error: the azide drawn with retention after an SN2 displacement
    route["compounds"]["7"]["smiles"] = "COC(=O)[C@@H]1C[C@H](N=[N+]=[N-])CN1C(=O)OC(C)(C)C"
    return route


def demo_items(language: str = "en") -> list[dict[str, Any]]:
    text = _TEXT.get(language, _TEXT["en"])
    image = _source_image(home_dir() / "demo")
    second: dict[str, Any] = {
        "id": "2",
        "label": text["l2"],
        "smiles": MISREAD_PROLINE,
        "source": text["src"],
        "note": text["n2"],
        "question": text["q2"],
        "candidates": [{"label": text["alt"], "smiles": TRUE_PROLINE}],
    }
    if image:
        second["image_path"] = image
    return [
        {"id": "1", "label": text["l1"], "smiles": "COC(=O)[C@@H](C)CO", "source": text["src"]},
        second,
        {"id": "3", "label": text["l3"], "smiles": "CC(O)C(C)C=CC", "source": text["src"], "note": text["n3"]},
    ]


def review_demo(language: str = "en", detached: bool = True, timeout_s: int = 7200, which: str = "both") -> dict[str, Any]:
    text = _TEXT.get(language, _TEXT["en"])
    sessions: dict[str, Any] = {}
    if which in {"both", "structures"}:
        sessions["structures"] = open_review(
            items=demo_items(language), detached=detached, timeout_s=timeout_s,
            title=text["s_title"], intro=text["s_intro"], language=language,
        )
    if which in {"both", "route"}:
        route = demo_route()
        route["title"] = text["r_title"]
        sessions["route"] = open_review(
            route=route, detached=detached, timeout_s=timeout_s, intro=text["r_intro"], language=language
        )
    return {
        "status": "ok",
        "sessions": {
            name: {key: session.get(key) for key in ("session_id", "review_url", "browser_opened", "expires_at")}
            for name, session in sessions.items()
        },
        "next_action": (
            "Two practice pages are open. Collect each result with "
            "`codex-chem review-result --session-id <id> --wait`."
        ),
    }
