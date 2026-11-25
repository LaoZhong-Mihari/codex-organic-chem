from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_skill_description_fits_strict_frontmatter_parsers():
    # The description is what makes a host load the skill. Keep it within the
    # Agent Skills limit of 1024 characters and a valid plain YAML scalar: an
    # unquoted ": " or " #" breaks strict parsers even where Claude Code copes.
    text = (REPOSITORY_ROOT / ".claude/skills/organic-chemistry-assistant/SKILL.md").read_text()
    frontmatter = text.split("---")[1]
    line = next(row for row in frontmatter.splitlines() if row.startswith("description: "))
    value = line.removeprefix("description: ")
    assert 0 < len(value) <= 1024
    assert ": " not in value and " #" not in value


def test_codex_and_claude_skill_instructions_are_identical():
    codex_skill = REPOSITORY_ROOT / ".agents/skills/organic-chemistry-assistant"
    claude_skill = REPOSITORY_ROOT / ".claude/skills/organic-chemistry-assistant"
    codex_files = {
        p.relative_to(codex_skill)
        for p in [codex_skill / "SKILL.md", *codex_skill.glob("references/**/*"), *codex_skill.glob("assets/**/*")]
        if p.is_file()
    }
    claude_files = {
        p.relative_to(claude_skill)
        for p in [claude_skill / "SKILL.md", *claude_skill.glob("references/**/*"), *claude_skill.glob("assets/**/*")]
        if p.is_file()
    }
    assert codex_files == claude_files
    for rel in sorted(codex_files):
        assert (codex_skill / rel).read_bytes() == (claude_skill / rel).read_bytes()
