from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_codex_and_claude_skill_instructions_are_identical():
    codex_skill = REPOSITORY_ROOT / ".agents/skills/organic-chemistry-assistant"
    claude_skill = REPOSITORY_ROOT / ".claude/skills/organic-chemistry-assistant"
    files = [codex_skill / "SKILL.md", *codex_skill.glob("references/**/*"), *codex_skill.glob("assets/**/*")]
    for source in files:
        if source.is_file():
            assert source.read_bytes() == (claude_skill / source.relative_to(codex_skill)).read_bytes()
