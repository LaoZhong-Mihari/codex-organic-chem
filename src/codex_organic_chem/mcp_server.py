from __future__ import annotations

import sys
from typing import Any

from .service import (
    chem_compute as _chem_compute,
    chem_draw as _chem_draw,
    chem_figure_tool_status as _chem_figure_tool_status,
    chem_input_review as _chem_input_review,
    chem_literature_search as _chem_literature_search,
    chem_mechanism_draft as _chem_mechanism_draft,
    chem_mechanism_render as _chem_mechanism_render,
    chem_mechanism_spec_example as _chem_mechanism_spec_example,
    chem_normalize_structure as _chem_normalize_structure,
    chem_ocsr_benchmark as _chem_ocsr_benchmark,
    chem_parse_image as _chem_parse_image,
    chem_parse_scheme as _chem_parse_scheme,
    chem_reaction_analyze as _chem_reaction_analyze,
    chem_synthesis_suggest as _chem_synthesis_suggest,
    chem_tool_doctor as _chem_tool_doctor,
)


def main() -> None:
    try:
        from mcp.server.fastmcp import FastMCP
    except Exception as exc:  # pragma: no cover - optional dependency
        sys.stderr.write(
            "The MCP dependency is not installed. Run `uv sync --extra mcp` "
            f"or install `mcp`. Import error: {exc}\n"
        )
        raise SystemExit(1)

    mcp = FastMCP("codex-organic-chem")

    @mcp.tool()
    def chem_parse_image(path: str, kind: str = "auto") -> dict[str, Any]:
        return _chem_parse_image(path=path, kind=kind)

    @mcp.tool()
    def chem_ocsr_benchmark(gold_smiles: str, image_dir: str) -> dict[str, Any]:
        return _chem_ocsr_benchmark(gold_smiles=gold_smiles, image_dir=image_dir)

    @mcp.tool()
    def chem_parse_scheme(image: str, crops: str, gold_map: str | None = None) -> dict[str, Any]:
        return _chem_parse_scheme(image=image, crops=crops, gold_map=gold_map)

    @mcp.tool()
    def chem_input_review(
        smiles: str | None = None,
        reaction_smiles: str | None = None,
        molfile: str | None = None,
        image_path: str | None = None,
        kind: str = "auto",
    ) -> dict[str, Any]:
        return _chem_input_review(
            smiles=smiles,
            reaction_smiles=reaction_smiles,
            molfile=molfile,
            image_path=image_path,
            kind=kind,
        )

    @mcp.tool()
    def chem_normalize_structure(smiles: str | None = None, molfile: str | None = None) -> dict[str, Any]:
        return _chem_normalize_structure(smiles=smiles, molfile=molfile)

    @mcp.tool()
    def chem_draw(
        smiles: str | None = None,
        reaction_smiles: str | None = None,
        output: str = "svg",
        output_file: str | None = None,
    ) -> dict[str, Any]:
        return _chem_draw(smiles=smiles, reaction_smiles=reaction_smiles, output=output, output_file=output_file)

    @mcp.tool()
    def chem_compute(
        smiles: str,
        tasks: list[str] | None = None,
        num_confs: int = 8,
        max_iters: int = 200,
    ) -> dict[str, Any]:
        return _chem_compute(smiles=smiles, tasks=tasks, num_confs=num_confs, max_iters=max_iters)

    @mcp.tool()
    def chem_tool_doctor() -> dict[str, Any]:
        return _chem_tool_doctor()

    @mcp.tool()
    def chem_figure_tool_status() -> dict[str, Any]:
        return _chem_figure_tool_status()

    @mcp.tool()
    def chem_literature_search(query: str, rows: int = 5) -> dict[str, Any]:
        return _chem_literature_search(query=query, rows=rows)

    @mcp.tool()
    def chem_reaction_analyze(input: str, mode: str = "sanity_check") -> dict[str, Any]:
        return _chem_reaction_analyze(input=input, mode=mode)

    @mcp.tool()
    def chem_mechanism_draft(
        reaction: str,
        style: str = "stepwise",
        quality: str = "draft",
        structure_confirmed: bool = False,
    ) -> dict[str, Any]:
        return _chem_mechanism_draft(
            reaction=reaction,
            style=style,
            quality=quality,
            structure_confirmed=structure_confirmed,
        )

    @mcp.tool()
    def chem_mechanism_render(spec: dict[str, Any], output_dir: str | None = None) -> dict[str, Any]:
        return _chem_mechanism_render(spec=spec, output_dir=output_dir)

    @mcp.tool()
    def chem_mechanism_spec_example() -> dict[str, Any]:
        return _chem_mechanism_spec_example()

    @mcp.tool()
    def chem_synthesis_suggest(
        target_smiles: str,
        stage: str = "first_disconnection",
        selected_option: str | None = None,
        confirmed: bool = False,
        literature: bool = True,
        literature_rows: int = 4,
    ) -> dict[str, Any]:
        return _chem_synthesis_suggest(
            target_smiles=target_smiles,
            stage=stage,
            selected_option=selected_option,
            confirmed=confirmed,
            literature=literature,
            literature_rows=literature_rows,
        )

    mcp.run()


if __name__ == "__main__":
    main()
