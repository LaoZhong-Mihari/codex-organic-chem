from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_codex_and_claude_skill_instructions_are_identical():
    codex_skill = REPOSITORY_ROOT / ".agents/skills/organic-chemistry-assistant/SKILL.md"
    claude_skill = REPOSITORY_ROOT / ".claude/skills/organic-chemistry-assistant/SKILL.md"

    assert codex_skill.read_bytes() == claude_skill.read_bytes()
