# Skills and MCP server

[Back to README](../README.md)

## Skill Setup (Codex and Claude Code)

Two synchronized copies of the skill ship with the project:

```text
.agents/skills/organic-chemistry-assistant/   # Codex
├── SKILL.md                 # the loop and the claim-to-check table
├── references/              # route-planning, reaction-outcome, structure-input,
│                            # human-review, figures, kit-api, mechanism-spec,
│                            # mechanism-graphics
├── assets/                  # example mechanism specs
└── agents/openai.yaml
.claude/skills/organic-chemistry-assistant/   # Claude Code (identical copy)
.claude/agents/chem-route-critic.md           # optional independent reviewer
```

`SKILL.md` is deliberately short; the task protocols in `references/` are read
only when the task needs them. Its `description` is what makes the host load
the skill at all; the README's benchmark section shows how often that happens.
Keep one copy per host: an older copy installed elsewhere under the same name
(for example in the user's skills or a synced account skill) is listed next to
this one with its own description, and the model may load either.

Claude Code discovers the `.claude/skills` copy automatically when this
repository is the working directory, and `.mcp.json` exposes the
`codex-organic-chem` MCP server to it. For any other environment, copy or
symlink the relevant folder into the target skills directory and keep this
repository available so the MCP server can launch (the Codex configuration is
in [INSTALL.md](../INSTALL.md#4-connect-the-mcp-server-and-the-skill)):

```bash
uv run --directory /path/to/this/repository --extra mcp codex-chem-mcp
```

## MCP

The repository includes `.mcp.json`, which starts the server with:

```bash
uv run --extra mcp codex-chem-mcp
```

`--extra mcp` makes uv install the `mcp` package on first start when the
environment was synced without it.

On connection the server sends the client short instructions: what the tools
check, and to load the `organic-chemistry-assistant` skill, when the host has
it, before the first tool call. `CODEX_CHEM_MCP_INSTRUCTIONS=off` leaves them
out.

The server registers 26 tools, plus `chem_python` when
`CODEX_CHEM_ENABLE_PYTHON=1`:

| Purpose | Tools |
|---|---|
| Checks that can fail | `chem_compare_structures`, `chem_enumerate`, `chem_route_check`, `chem_model_kit` |
| Structures and properties | `chem_normalize_structure`, `chem_draw`, `chem_compute` (RDKit, GFN2-xTB, CREST) |
| Reading images | `chem_parse_image`, `chem_parse_scheme`, `chem_ocsr_benchmark` |
| Human review | `chem_input_review`, `chem_structure_review_batch`, `chem_structure_review_result`, `chem_route_review`, `chem_review_mode` |
| Reactions and routes | `chem_reaction_analyze`, `chem_synthesis_suggest`, `chem_literature_search` |
| Mechanisms | `chem_mechanism_draft`, `chem_mechanism_validate`, `chem_mechanism_render`, `chem_mechanism_spec_example` |
| Figures | `chem_route_figure`, `chem_route_figure_spec_example`, `chem_figure_tool_status` |
| Set-up | `chem_tool_doctor` |
| Opt-in | `chem_python` (Python with RDKit in a temporary directory, for hosts without a shell) |
