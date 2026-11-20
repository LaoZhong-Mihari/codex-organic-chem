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
    files = [codex_skill / "SKILL.md", *codex_skill.glob("references/**/*"), *codex_skill.glob("assets/**/*")]
    for source in files:
        if source.is_file():
            assert source.read_bytes() == (claude_skill / source.relative_to(codex_skill)).read_bytes()
